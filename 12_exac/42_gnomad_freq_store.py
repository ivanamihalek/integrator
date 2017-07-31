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
			print field
			for k in storable_fields:
			#for k in named_csqs.keys():
				print k, ":",  named_csqs[k]
				print "-------------"
		# somebody screwed up in gnoMAD so the leading base is not included in insert or deletion
		# what if there is some more complex variant?
		if False and  len(ref)>1 and len(named_csqs['Allele'])>1 and abs(len(ref)-len(named_csqs['Allele']))>1:
			csq =  named_csqs['Consequence']
			if csq in ['insertion', 'deletion'] and  not 'non_coding' in csq and ('splice' in csq  or 'exon' in csq or 'coding'  in csq):
				# this is only worth so much to me at this point; mental note: do your own annotation
				probably_allele = ref[0]+named_csqs['Allele']
				allele_matches = [var for var in variants if var==probably_allele]
				if len(allele_matches) != 1:
					print "point   {}  {}  {} ---> {}   {}".format(named_csqs['VARIANT_CLASS'], named_csqs['Consequence'], ref, named_csqs['Allele'], vars)
					#exit()
				# the annotation is complete garbage
				# I see thins like  deletion    GC ---> T   G,TC,AC
				# it does seem to be limited to up/downstream variants; I have not found exons or splice sites suffering from this
				# so as some sort of a patch, if the allele with prepended leading base is not found among the variants, leave the csq field empty
		patched_allele = None

		if named_csqs['VARIANT_CLASS']=='SNV':
			# this is ok, G-->T and similar
			# somtimes I just don't know hot to reconcile this:
			# Example: ref variant  CT	C,GT
			# and the nnotation says
			# G|intron_variant|MODIFIER|RPS4Y1|ENSG00000129824|Transcript|ENST00000250784|protein_coding||1/6|ENST00000250784.8:c.3+17C>G||||||||2||1||deletion
			# deletion ?!
			# there is HGNC recommended annotation - look a the whole thing at some other time;
			# ideally just throw the whole crap away and do my own annotation
			patched_allele = named_csqs['Allele'][0] + ref[1:]
		elif named_csqs['VARIANT_CLASS']=='insertion':
			patched_allele = ref[0]+named_csqs['Allele']
		elif named_csqs['VARIANT_CLASS']=='deletion':
			patched_allele = ref[0]
		# exactly one match in variants
		if len([var for var in variants if var==patched_allele])==1:
			named_csqs['Allele'] = patched_allele
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
			for g in gender+['both']: data[cnt][p][g] = 0
	consequences = ""
	for  k,v  in named_fields.iteritems():
		field = k.split("_")
		if field[0]=='CSQ':
			consequences = parse_csq(v, consequence_header_fields, ref, vars)
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
	consequences, counts = parse_info (chrom, ref, variants, info_string, consequence_header_fields)
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
def main():
	infile = open("/databases/exac/gnomad.exomes.r2.0.1.sites.vcf")
	#infile = open("/databases/exac/gnomad_test.txt")
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
		chrom = line.split("\t")[0]
		#if not chrom in ['1','9','10','22']: continue
		#if not chrom in ['2','8','11','21']: continue
		#if not chrom in ['3','7','12','20']: continue
		#if not chrom in  ['4','6', '13','19']: continue
		#if not chrom in  ['5','14','16','18']: continue
		if not chrom in  ['15','17', 'X','Y']: continue
		ret = process_line(line, consequence_header_fields)
		count += 1
		#if count%100000==0: print "%dK lines out of 15014K (%5.2f%%)" % (count/1000,  float(count)/15014744*100)
		if count%10000==0: print "%dK lines of  %s" % (count/1000,  chrom)
		if not ret: continue
		[chrom, addr, ref, variants, consequences,  ac, an, ac_afr, an_afr,
		 ac_amr, an_amr, ac_asj, an_asj, ac_eas, an_eas, ac_fin, an_fin,
		 ac_nfe, an_nfe, ac_oth, an_oth, ac_sas, an_sas] = ret
		if False and (len(ref)>1 or len(variants)>1):
			print ref, variants # <<<<<<<<<<<<<<<<  !!!!!!!!!!!!
			print chrom, addr
			print "\n".join(consequences.split(","))
			print

		table = "gnomad_freqs_chr_" + chrom
		fixed_fields  = {'position':int(addr)}
		update_fields = {'reference':ref, 'variants':variants, 'variant_counts':ac, 'total_count':an}
		update_fields['consequences'] = consequences
		update_fields['afr_counts'] = ac_afr
		update_fields['afr_tot_count'] = an_afr
		update_fields['amr_counts'] = ac_amr
		update_fields['amr_tot_count'] = an_amr
		update_fields['asj_counts'] = ac_asj
		update_fields['asj_tot_count'] = an_asj
		update_fields['eas_counts'] = ac_eas
		update_fields['eas_tot_count'] = an_eas
		update_fields['fin_counts'] = ac_fin
		update_fields['fin_tot_count'] = an_fin
		update_fields['nfe_counts'] = ac_nfe
		update_fields['nfe_tot_count'] = an_nfe
		update_fields['oth_counts'] = ac_oth
		update_fields['oth_tot_count'] = an_oth
		update_fields['sas_counts'] = ac_sas
		update_fields['sas_tot_count'] = an_sas
		store_or_update(cursor, table, fixed_fields, update_fields, verbose=False, primary_key='position')
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914