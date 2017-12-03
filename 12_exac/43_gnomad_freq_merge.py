#!/usr/bin/python

from integrator_utils.mysql import *
import shlex

count_fields = ["AC", "AN", "AF", "GC", "Hemi"]
population_fields = ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE","OTH","SAS"]
gender = ["Male","Female"]
# lof fields seem to be mostly empty
storable_fields = ['Allele', 'SYMBOL', 'SWISSPROT', 'VARIANT_CLASS', 'Consequence', 'cDNA_position', 'CDS_position',
				   'Protein_position', 'Amino_acids','SIFT','PolyPhen'] # 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']
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
					print "point   {}  {}  {} ---> {}   {}".format(named_csqs['VARIANT_CLASS'], named_csqs['Consequence'], ref, named_csqs['Allele'], vars)
					#exit()
				# the annotation is complete garbage
				# I see things like  deletion    GC ---> T   G,TC,AC
				# it does seem to be limited to up/downstream variants; I have not found exons or splice sites suffering from this
				# so as some sort of a patch, if the allele with prepended leading base is not found among the variants, leave the csq field empty
		patched_allele = None

		if named_csqs['Allele']=='-':
			patched_allele = ref[0]
		elif len(ref)==1 and len(named_csqs['Allele'])==1:
			patched_allele = named_csqs['Allele']
		elif 'insertion' in csq or named_csqs['VARIANT_CLASS']=='insertion':
			patched_allele = ref[0]+named_csqs['Allele']
		elif 'deletion' in csq or named_csqs['VARIANT_CLASS']=='deletion':
			patched_allele = ref[0]
		elif 'missense' in csq or 'stop' in csq: # this should presumably be an SNV; there are deletions marked as SNV
			# this is ok, G-->T and similar
			# somtimes I just don't know hot to reconcile this:
			# Example: ref variant  CT	C,GT
			# and the nnotation says
			# G|intron_variant|MODIFIER|RPS4Y1|ENSG00000129824|Transcript|ENST00000250784|protein_coding||1/6|ENST00000250784.8:c.3+17C>G||||||||2||1||deletion
			# deletion ?!
			# there is HGNC recommended annotation - look a the whole thing at some other time;
			# ideally just throw the whole crap away and do my own annotation
			patched_allele = named_csqs['Allele'][0] + ref[1:]

		# exactly one match in variants
		if patched_allele and  len([var for var in variants if var==patched_allele])==1:
			named_csqs['Allele'] = patched_allele
			csq = "|".join([named_csqs[k] for k in storable_fields])
			if not csq in consequences: consequences.append(csq)

	#print "***************************"
	#print ",".join(consequences)
	#print "***************************"
	#print

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
	for k,v in named_fields.iteritems():
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
def process_line(line, consequence_header_fields):
	fields = line.rstrip().split("\t")
	[chrom, addr, junk, ref, variants, quality_score,  filter, info_string] = fields[:8]
	#if chrom != '3':  return None # <<<<<<<<<<<<<<<<  !!!!!!!!!!!!
	if filter!='PASS': return None
	consequences, counts = parse_info(chrom, ref, variants, info_string, consequence_header_fields)
	if not consequences: return None
	return [chrom, addr, ref, variants, consequences] + counts


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
def merge(cursor, table, col_names, new_fields, verbose):
	qry  = "select * from %s where position=%d" % (table, new_fields['position'])
	rows = search_db(cursor, qry, verbose)
	old_fields =  dict(zip(col_names, rows[0])) # there should be only one position, with all variants listed here
	# note: the reference may be different!

	if new_fields['variants']==old_fields['variants']: return
	new_vars_nominal = new_fields['variants'].split(',')
	old_vars = old_fields['variants'].split(',')
	if len(old_vars)<3 and len(new_vars_nominal)<3: return
	new_vars = [str(v) for v in new_vars_nominal if not v in old_vars]
	if len(new_vars)==0: return
	print "merging"
	print "new_vars_nominal", new_vars_nominal
	print "old_vars", old_vars
	print "new_vars", new_vars
	print "new_fields['reference']", new_fields['reference']
	print "consequences", new_fields['consequences']
	for v in new_vars:
		i = new_vars_nominal.index(v)
		print "i, v ", i, v
		old_vars.append(v)
		# parse counts
		for f in ['variant_counts', 'afr_counts', 'afr_tot_count', 'amr_counts', 'amr_tot_count',
		          'asj_counts', 'asj_tot_count', 'eas_counts',  'eas_tot_count', 'fin_counts',
		          'fin_tot_count', 'nfe_counts', 'nfe_tot_count',
		          'oth_counts', 'oth_tot_count', 'sas_counts', 'sas_tot_count']:
			if ',' in  new_fields[f]:
				new_val = str(new_fields[f]).split(',')[i]
			elif new_fields[f]=='0':
				new_val = '0'
			else:
				new_val = new_fields[f]
			print f, new_val

		# parse consequences
		new_cons_list = []
		for cons in  new_fields['consequences'].split(","):
			if cons.split('|')[0].replace(" ","") == v:
				new_cons_list.append(cons)
		if len(new_cons_list)>0:
			new_cons = new_cons_list.join(",")
		else:
			new_cons = None

		print 'consequences', new_cons

	exit()

	return


##########################################
def store_or_merge(cursor, table, col_names, fields, verbose):

	# check if the row exists
	qry  = "select position from %s where position=%d"  % (table, fields['position'])
	rows = search_db(cursor, qry, verbose=False)

	exists = rows and (type(rows[0][0]) is long)
	if exists:
		merge(cursor, table, col_names, fields, verbose)
	else:
		#print all_fields
		#store_without_checking(cursor, table, fields, verbose)
		pass
##########################################
def get_col_names(cursor, table):
	qry  = "SELECT `COLUMN_NAME`  FROM `INFORMATION_SCHEMA`.`COLUMNS` "
	qry += "WHERE `TABLE_SCHEMA`='blimps_development' "
	qry += "AND `TABLE_NAME`='%s' " % table
	rows = search_db(cursor, qry, verbose=True)
	col_names = [row[0] for row in rows]
	return col_names

##########################################
def main():
	#infile = open("/databases/exac/gnomad_test.txt")

	db, cursor = connect()

	reading = False
	count = 0
	consequence_header_fields = None

	chrom = "15"
	infile = open("/databases/exac/release_2.0.2_vcf_genomes_gnomad.genomes.r2.0.2.sites.chr%s.vcf" % chrom)
	table = "gnomad_freqs_chr_" + chrom
	col_names = get_col_names(cursor, table)

	for line in infile:
		if not reading:
			if line[:6] == '#CHROM': reading=True
			if not consequence_header_fields  and line[:len('##INFO=<')] == '##INFO=<':
				consequence_header_fields = find_csq_header_fields(line)
			continue
		if not consequence_header_fields:
			print 'CSQ info not found'
			exit()
		chrom_inline = line.split("\t")[0]
		if chrom_inline!=chrom:
			print "chromosome %s?" % chrom
			print line
			exit()
		ret = process_line(line, consequence_header_fields)

		count += 1
		#if count%100000==0: print "%dK lines out of 15014K (%5.2f%%)" % (count/1000,  float(count)/15014744*100)
		if count%10000==0: print "%dK lines of  %s" % (count/1000,  chrom)
		if not ret: continue # error or silent mutation
		[chrom, addr, ref, variants, consequences,  ac, an, ac_afr, an_afr,
		 ac_amr, an_amr, ac_asj, an_asj, ac_eas, an_eas, ac_fin, an_fin,
		 ac_nfe, an_nfe, ac_oth, an_oth, ac_sas, an_sas] = ret


		if False and (len(ref)>1 or len(variants)>1):
			print ref, variants # <<<<<<<<<<<<<<<<  !!!!!!!!!!!!
			print chrom, addr
			print "\n".join(consequences.split(","))
			print

		all_fields = {'position':int(addr), 'reference':ref, 'variants':variants, 'variant_counts':ac, 'total_count':an}

		all_fields['consequences'] = consequences
		all_fields['afr_counts'] = ac_afr
		all_fields['afr_tot_count'] = an_afr
		all_fields['amr_counts'] = ac_amr
		all_fields['amr_tot_count'] = an_amr
		all_fields['asj_counts'] = ac_asj
		all_fields['asj_tot_count'] = an_asj
		all_fields['eas_counts'] = ac_eas
		all_fields['eas_tot_count'] = an_eas
		all_fields['fin_counts'] = ac_fin
		all_fields['fin_tot_count'] = an_fin
		all_fields['nfe_counts'] = ac_nfe
		all_fields['nfe_tot_count'] = an_nfe
		all_fields['oth_counts'] = ac_oth
		all_fields['oth_tot_count'] = an_oth
		all_fields['sas_counts'] = ac_sas
		all_fields['sas_tot_count'] = an_sas
		store_or_merge(cursor, table, col_names, all_fields, verbose=False)
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914