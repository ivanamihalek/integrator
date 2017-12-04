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
	print fields
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
	print ">>>>>>>>>>>"
	print named_csqs['Consequence']
	print named_csqs['Allele']
	print "<<<<<<<<<<<"
	print "consequences1", consequences
	return ",".join(consequences)

#########################################
def parse_info(chrom, ref, vars, info, consequence_header_fields):

	fields = info.split(";")
	print " ***", fields
	named_fields= dict (map(lambda x: x.split("="), [f for f in fields if '=' in f]))
	data = {}
	for cnt in count_fields:
		data[cnt] = {}
		for p in population_fields+['overall']:
			data[cnt][p] = {}
			for g in gender+['both']: data[cnt][p][g] = 0
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

	print "consequences2", consequences
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
	print info_string
	print
	print filter
	print
	if filter!='PASS': return None
	consequences, counts = parse_info (chrom, ref, variants, info_string, consequence_header_fields)
	print  "consequences3", consequences
	if not consequences: return None
	return [chrom, addr, ref, variants, consequences] + counts

##########################################
def store_variant(cursor, table, var_fields_unpacked):

	print "var_fields_unpacked", var_fields_unpacked
	print
	[addr, ref, variant, consequences,  ac, an, ac_afr, an_afr,
	 ac_amr, an_amr, ac_asj, an_asj, ac_eas, an_eas, ac_fin, an_fin,
	 ac_nfe, an_nfe, ac_oth, an_oth, ac_sas, an_sas] = var_fields_unpacked
	
	fixed_fields  = {'position':int(addr)}
	
	update_fields = {'reference':ref, 'variant':variant, 'variant_count':int(ac), 'total_count':int(an)}
	update_fields['consequences'] = consequences
	update_fields['afr_count'] = int(ac_afr)
	update_fields['afr_tot_count'] = int(an_afr)
	update_fields['amr_count'] = int(ac_amr)
	update_fields['amr_tot_count'] = int(an_amr)
	update_fields['asj_count'] = int(ac_asj)
	update_fields['asj_tot_count'] = int(an_asj)
	update_fields['eas_count'] = int(ac_eas)
	update_fields['eas_tot_count'] = int(an_eas)
	update_fields['fin_count'] = int(ac_fin)
	update_fields['fin_tot_count'] = int(an_fin)
	update_fields['nfe_count'] = int(ac_nfe)
	update_fields['nfe_tot_count'] = int(an_nfe)
	update_fields['oth_count'] = int(ac_oth)
	update_fields['oth_tot_count'] = int(an_oth)
	update_fields['sas_count'] = int(ac_sas)
	update_fields['sas_tot_count'] = int(an_sas)
	store_or_update(cursor, table, fixed_fields, update_fields, verbose=False)

	return

##########################################
def process_line(cursor, line, consequence_header_fields):

	ret = parse_data_line(line, consequence_header_fields)

	if not ret: return # error or silent mutation
	keys = ['chrom', 'addr', 'ref', 'variants', 'consequences',  'ac', 'an', 'ac_afr', 'an_afr',
	 'ac_amr', 'an_amr', 'ac_asj', 'an_asj', 'ac_eas', 'an_eas', 'ac_fin', 'an_fin',
	 'ac_nfe', 'an_nfe', 'ac_oth', 'an_oth', 'ac_sas', 'an_sas']
	named_fields =  dict(zip(keys,ret))


	# unpack variants

	split_string = {}
	split_string['variants'] = named_fields['variants'].split(',')
	number_of_variants = len(split_string['variants'])

	for k in  ['ac', 'an', 'ac_afr', 'an_afr',
	 'ac_amr', 'an_amr', 'ac_asj', 'an_asj', 'ac_eas', 'an_eas', 'ac_fin', 'an_fin',
	 'ac_nfe', 'an_nfe', 'ac_oth', 'an_oth', 'ac_sas', 'an_sas']:
		val = named_fields[k]
		split_string[k] = val.split(',')
		if len(split_string[k])!=number_of_variants:
			if named_fields[k]=='0':
				print 'err <<< '
				exit()
			else:
				split_string[k] = ['0']*number_of_variants

	table = "gnomad_freqs_chr_" + named_fields['chrom']
	
	print split_string['variants']
	for i in range(len(split_string['variants'])):
		consq_for_this_variant = []
		v = split_string['variants'][i]
		for c in named_fields['consequences'].split(','):
			if c.split('|')[0] == v: consq_for_this_variant.append(c)
			var_fields_unpacked = [named_fields['addr'],
			    named_fields['ref'], v, ",".join(consq_for_this_variant)]
			for k in  ['ac', 'an', 'ac_afr', 'an_afr',
				 'ac_amr', 'an_amr', 'ac_asj', 'an_asj', 'ac_eas', 'an_eas', 'ac_fin', 'an_fin',
				 'ac_nfe', 'an_nfe', 'ac_oth', 'an_oth', 'ac_sas', 'an_sas']:
				var_fields_unpacked.append(split_string[k][i])
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
	#infile = open("/databases/exac/release_2.0.2_vcf_exomes_gnomad.exomes.r2.0.2.sites.vcf")
	infile = open("/databases/exac/gnomad_test.txt")
	#infile = open("/databases/exac/test_pccb.vcf") # <<<<<<<<<<<<<<<<  !!!!!!!!!!!!
	#infile = open("/databases/exac/testY.vcf")

	db, cursor = connect()

	reading = False
	count = 0
	consequence_header_fields = None
	chrom = ""
	for line in infile:
		if not reading:
			if line[:6] == '#CHROM': reading=True
			if not consequence_header_fields  and line[:len('##INFO=<')] == '##INFO=<':
				consequence_header_fields = find_csq_header_fields(line)
			continue
		if not consequence_header_fields:
			print 'CSQ info not found'
			exit()
		count += 1
		if count%10000==0: print "%dK lines of  %s" % (count/1000,  chrom)

		chrom = line.split("\t")[0]
		#if not chrom in ['1','9','10','22']: continue
		#if not chrom in ['2','8','11','21']: continue
		#if not chrom in ['3','7','12','20']: continue
		#if not chrom in  ['4','6', '13','19']: continue
		#if not chrom in  ['5','14','16','18']: continue
		#if not chrom in  ['15','17', 'X','Y']: continue
		process_line(cursor,line, consequence_header_fields)

	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()

# 8,608,914