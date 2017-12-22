#!/usr/bin/python
# create index gnomad_pos_idx_1 on gnomad_freqs_chr_1 (position);
# careful with this, because it adds to existing without checking
# perhaps some better wasy could be found for this
from integrator_utils.mysql import *
import shlex

count_fields = ["AC", "AN", "AF", "GC", "Hemi"]
population_fields = ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE","OTH","SAS"]
gender = ["Male","Female"]
# lof fields seem to be mostly empty
storable_fields = ['Allele', 'SYMBOL', 'SWISSPROT', 'VARIANT_CLASS', 'Consequence', 'cDNA_position', 'CDS_position',
                   'Protein_position', 'Amino_acids','SIFT','PolyPhen'] # 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']

all_keys = ['chrom', 'addr', 'ref', 'variants', 'consequences',  'ac', 'an', 'ac_afr', 'an_afr',
	 'ac_amr', 'an_amr', 'ac_asj', 'an_asj', 'ac_eas', 'an_eas', 'ac_fin', 'an_fin',
	 'ac_nfe', 'an_nfe', 'ac_oth', 'an_oth', 'ac_sas', 'an_sas']
counts = [key for key in all_keys if key[:2]=='ac']
totals = [key for key in all_keys if key[:2]=='an']

db_gnomad_freq_keys = []

# the annotation is crap ...
# 22:17669352 G>GGGACA is an SNV in somebody's book


#########################################
def check_difference (ref,alt): # bcs of the way the info is store, can have things like "ATT, CTT", inidcating A->C
	# this is untested
	if ref==alt: return ref,alt
	i=-1
	while True:
		if not ref[i]: break
		if not alt[i]: break
		if ref[i] != alt[i]: break
		i -= 1

	frm = ref[0:i+1]
	to  = alt[0:i+1]
	if not frm or frm=="": frm = "-"
	if not to  or to=="": to = "-"
	return frm,to

#########################################
def parse_csq(csq_string, consequence_header_fields, ref, vars):
	fields = csq_string.split(",")
	variants = vars.split(",")
	consequences = []
	for field in fields:
		named_csqs = dict(zip(consequence_header_fields,field.split("|")))
		#if 'missense' in named_csqs['Consequence']:
		if False:
			print
			print "ref:",ref
			print "annotation filed:", field
			for k in storable_fields:
			#for k in named_csqs.keys():
				print k, ":",  named_csqs[k]
				print "-------------"

		if named_csqs['Consequence'] == 'intergenic_variant': return None

		csq =  named_csqs['Consequence']
		if 'synonymous' in csq:
			return None
		# somebody screwed up in gnoMAD so the leading base is not included in insert or deletion
		# what if there is some more complex variant?
		if False and  len(ref)>1 and len(named_csqs['Allele'])>1 and abs(len(ref)-len(named_csqs['Allele']))>1:

			if csq in ['insertion', 'deletion'] and  not 'non_coding' in csq and ('splice' in csq  or 'exon' in csq or 'coding'  in csq):
				# this is only worth so much to me at this point; mental note: do your own annotation
				probably_allele = ref[0]+named_csqs['Allele']
				allele_matches = [var for var in variants if var==probably_allele]
				if len(allele_matches) != 1:
					print "point   {}  {}  {} ---> {}  {}".format(named_csqs['VARIANT_CLASS'], named_csqs['Consequence'], ref, named_csqs['Allele'], vars)
					#exit()
				# the annotation does not follow the same variant representation as the variants field
		patched_allele = None

		# exactly one match in variants
		csq = "|".join([named_csqs[k] for k in storable_fields])
		if not csq in consequences: consequences.append(csq)

	return ",".join(consequences)

#########################################
def parse_info(chrom, ref, vars, info, consequence_header_fields):

	fields = info.split(";")
	named_fields= dict (map(lambda x: x.split("="), [f for f in fields if '=' in f]))
	data = {}
	for cnt in count_fields:
		data[cnt] = {}
		for p in population_fields+['overall']:
			data[cnt][p] = {}
			for g in gender+['both']: data[cnt][p][g] = '0'
	consequences = ""
	for  k,v  in named_fields.iteritems():
		field = k.split("_")
		if field[0]=='CSQ':
			consequences = parse_csq(v, consequence_header_fields, ref, vars)
			if not consequences: return None, None # synonymous mutation
			continue
		if not field[0] in count_fields: continue
		count_type = field[0]
		l = len(field)
		if l==1:
			data[count_type]['overall']['both'] = v
		elif l==2:
			if field[1] in gender:
				data[count_type]['overall'][field[1]] = v
			elif field[1] in population_fields:
				data[count_type][field[1]]['both'] = v
			elif field[1] in ['raw','POPMAX']:
				continue
			else:
				print "unrecognized key:", k
				exit()
		elif l==3:
			[p,g] = field[1:3]
			if not p in population_fields or not g in gender:
				print "unrecognized key:", k
				exit()
			data[count_type][p][g] = v
		else:
			print "unrecognized key:", k
			exit()

	# sanity check
	if False:
		if chrom in ['X','Y']:
			for p in population_fields:
				for g in gender+['both']:
					print p, g, data['AC'][p][g], data['AN'][p][g], data['AF'][p][g]
			print 'overall', data['AC']['overall']['both'], data['AN']['overall']['both'], data['AF']['overall']['both']
			print '-----'
		else:
			# AN = total number of alleles
			total = 0
			for p in population_fields:
				total += int(data['AN'][p]['both'])
				print p, data['AC'][p]['both'], data['AN'][p]['both'], data['AF'][p]['both']
			print data['AC']['overall']['both'], data['AN']['overall']['both'], data['AF']['overall']['both']
			if int(data['AN']['overall']['both']) != total:
				print "count mismatch:",data['AN']['overall']['both'], total
				exit()

	counts = []
	counts.append(data['AC']['overall']['both'])
	counts.append(data['AN']['overall']['both'])
	for p in population_fields:
		counts.append(data['AC'][p]['both'])
		counts.append(data['AN'][p]['both'])
	return consequences, counts

#########################################
def parse_data_line(line, consequence_header_fields):
	fields = line.rstrip().split("\t")
	[chrom, addr, junk, ref, variants, quality_score,  filter, info_string] = fields[:8]
	if filter!='PASS': return None
	consequences, counts = parse_info (chrom, ref, variants, info_string, consequence_header_fields)

	if not consequences: return None
	return [chrom, addr, ref, variants, consequences] + counts

##########################################
def minimal_var_representation(str1,str2):

	maxl = max([len(str1),len(str2)])
	# pad to the same length
	str1 += "x"*(maxl-len(str1))
	str2 += "x"*(maxl-len(str2))
	min1 = str1[:1]
	min2 = str2[:1]
	for i in range(1,maxl):
		if str1[i]==str2[i]: break
		min1 += str1[i]
		min2 += str2[i]

	return min1.replace("x",""), min2.replace("x","")

##########################################
def store_or_add_variant(cursor, table, fixed_fields, update_fields, verbose=False):

	qry  = "select * from %s "   % table
	qry += "where position=%d "  % fixed_fields['position']
	qry += "and reference='%s' " % fixed_fields['reference']
	qry += "and variant='%s' "   % fixed_fields['variant']
	ret = search_db(cursor,qry)
	if ret and ret[0][1]==fixed_fields['position']:
		existing_fields = dict(zip(db_gnomad_freq_keys,ret[0]))
		# sum and update
		update_assignments = []
		# the first 5 are addres, consequences etc, the last one is hotspot id
		for count in db_gnomad_freq_keys[5:-1]:
			#print count, existing_fields[count], update_fields[count]
			update_fields[count] += existing_fields[count]
			# create update query
			update_assignments.append("%s=%d" %(count,update_fields[count]))
		# create update query
		qry = "update %s set " % table
		qry += ", ".join(update_assignments) + " "
		qry += "where id=%d "  % existing_fields['id']
		# print update qry
		# print qry
		# print fixed_fields['position']
		# print fixed_fields['reference']
		# print fixed_fields['variant']
		search_db(cursor,qry)

	else:
		all_fields = fixed_fields
		all_fields.update(update_fields)
		# print "new entry"
		store_without_checking(cursor, table, all_fields, verbose=False)


##########################################
def store_variant(cursor, table, var_fields_unpacked):

	# not sure why this happens - filtered varaints (the ones that do not have 'PASS')
	# still may have count > 0
	# for the gnomad blogpost
	# "All sites with no high quality genotype (AC = 0) are marked as filtered using the AC0 filter."
	if var_fields_unpacked['ac']=='0':
		return

	# in trying for a compact storage, gnomad packs several variants into a single line
	# for example ref CT, variants C,AT to account for a deletion and an SNV
	# we are going to tun that back int CT, C  and C, A
	ref,var = minimal_var_representation(var_fields_unpacked['ref'], var_fields_unpacked['variant'])

	consq_for_this_variant = []
	# however in the consequences field, gnomad uses the minimalist representation
	if len(ref)>1 and len(var)==1:
		var_rep = "-"
	elif  len(ref)==1 and len(var)>1:
		var_rep = var[1:]
	else:
		var_rep = var
	for c in var_fields_unpacked['consequences'].split(','):
		csq_list = c.split('|')
		if csq_list[0] == var_rep: consq_for_this_variant.append('|'.join([var]+csq_list[1:]))

	fixed_fields  = {'position':int(var_fields_unpacked['addr']),'reference': ref,  'variant': var}
	update_fields = {}

	update_fields['consequences'] = ",".join(consq_for_this_variant)
	for count in counts:
		if count == 'ac':
			new_key = 'variant_count'
		else:
			population = count.split("_")[-1]
			new_key = "_".join([population,'count'])
		update_fields[new_key] = int(var_fields_unpacked[count])

	for count in totals:
		if count == 'an':
			new_key = 'total_count'
		else:
			population = count.split("_")[-1]
			new_key = "_".join([population,'tot','count'])
		update_fields[new_key] = int(var_fields_unpacked[count])

	store_or_add_variant(cursor, table, fixed_fields, update_fields, verbose=False)

	return

##########################################
def process_line(cursor, line, consequence_header_fields):

	ret = parse_data_line(line, consequence_header_fields)

	if not ret: return # error or silent mutation

	# unpack variants
	named_fields = dict(zip(all_keys,ret))
	split_string = {}
	split_string['variants'] = named_fields['variants'].split(',')
	number_of_variants = len(split_string['variants'])

	for k in  counts:
		val = named_fields[k]
		split_string[k] = val.split(',')
		if len(split_string[k])!=number_of_variants:
			if named_fields[k]!='0':
				print 'err <<< ', k, number_of_variants, split_string[k]
				exit()
			else:
				split_string[k] = ['0']*number_of_variants

	table = "gnomad_freqs_chr_" + named_fields['chrom']
	
	for i in range(len(split_string['variants'])):
		v = split_string['variants'][i]
		# consequences will be filtered fro the varaint they refer to in store_variant,
		# after we determine the minimal representation for the variant
		var_fields_unpacked = {'addr': named_fields['addr'],
			                    'ref': named_fields['ref'], 'variant': v,
			                    'consequences': named_fields['consequences']}
		for k in counts:
			var_fields_unpacked[k] = split_string[k][i]
		for k in totals:
			var_fields_unpacked[k] = named_fields[k]

		store_variant(cursor, table, var_fields_unpacked)
	return

#########################################
def find_csq_header_fields(line):
	my_splitter = shlex.shlex(line.rstrip()[len('##INFO=<'):-1], posix=True)
	my_splitter.whitespace += ','
	my_splitter.whitespace_split = True
	fields = list(my_splitter)
	field_name = fields[0].split('=')[-1]
	field_content = fields[-1].split('=')[-1]
	if field_name == 'CSQ':
		return field_content.split(':')[1].replace(" ","").split('|')
	return None

##########################################
def main():
	#infile = open("/databases/exac/gnomad.genomes.r2.0.2.sites.coding_only.chrX.vcf")
	infile = open("/databases/exac/gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.vcf")
	#infile = open("/databases/exac/gnomad_test.2.txt")
	#infile = open("/databases/exac/test_pccb.vcf") # <<<<<<<<<<<<<<<<  !!!!!!!!!!!!
	#infile = open("/databases/exac/testY.vcf")

	# This one is good. Form gnomad blogpost, literally:
	# "There is no chromosome Y coverage / calls in the gnomAD genomes, because reasons."

	db, cursor = connect()
	qry  = "SELECT `COLUMN_NAME`  FROM `INFORMATION_SCHEMA`.`COLUMNS` "
	qry += "WHERE `TABLE_SCHEMA`='blimps_development' AND `TABLE_NAME`='gnomad_freqs_chr_1'"

	global db_gnomad_freq_keys
	db_gnomad_freq_keys = [entry[0] for entry in search_db(cursor, qry)]


	reading = False
	count = 0
	consequence_header_fields = None
	chrom = ""

	# redo 22
	#chromgroup = ['1','9','10','22']
	#chromgroup  = ['2','8','11','21']
	#chromgroup  = ['3','7','12','22']
	#chromgroup = ['4','6', '13','19']
	#chromgroup = ['5','14','16','18']
	chromgroup = ['15','17', 'X','Y']
	chromgroup = ['20']
	chromgroup = ['22']

	print "chomosomes:", chromgroup

	for line in infile:
		if not reading:
			if line[:6] == '#CHROM': reading=True
			if not consequence_header_fields and line[:len('##INFO=<')] == '##INFO=<':
				consequence_header_fields = find_csq_header_fields(line)
			continue
		if not consequence_header_fields:
			print 'CSQ info not found'
			exit()

		chrom = line.split("\t")[0]
		if not chrom in chromgroup: continue
		count += 1
		if count%1000==0: print "%dK lines of  %s" % (count/1000,  chrom)
		process_line(cursor,line, consequence_header_fields)

	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()

# 8,608,914