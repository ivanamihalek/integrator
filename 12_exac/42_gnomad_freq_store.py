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
def parse_csq(csq_string, consequence_header_fields, ref):
	fields = csq_string.split(",")
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
		if len(ref)>1 and len(named_csqs['Allele'])>1:
				print "{} ---> {}".format(ref, named_csqs['Allele'])
		if named_csqs['VARIANT_CLASS']=='SNV':
			pass # this is ok, G-->T and similar
		elif named_csqs['VARIANT_CLASS']=='insertion':
			named_csqs['Allele'] = ref[0]+named_csqs['Allele']
		elif named_csqs['VARIANT_CLASS']=='deletion':
			named_csqs['Allele'] = ref[0]
		csq = "|".join([named_csqs[k] for k in storable_fields])
		if not csq in consequences: consequences.append(csq)
	return ",".join(consequences)

#########################################
def parse_info(chrom, ref, info, consequence_header_fields):

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
			consequences = parse_csq(v, consequence_header_fields, ref)
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
	consequences, counts = parse_info (chrom, ref, info_string, consequence_header_fields)
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
	#infile = open("/databases/exac/gnomad.exomes.r2.0.1.sites.vcf")
	#infile = open("/databases/exac/gnomad_test.txt")
	infile = open("/databases/exac/test_pccb.vcf") # <<<<<<<<<<<<<<<<  !!!!!!!!!!!!

	db, cursor = connect()

	reading = False
	count = 0
	consequence_header_fields = None
	for line in infile:
		if not reading:
			if line[:6] == '#CHROM': reading=True
			if not consequence_header_fields  and line[:len('##INFO=<')] == '##INFO=<':
				consequence_header_fields = find_csq_header_fields(line)
			continue
		if not consequence_header_fields:
			print 'CSQ info not found'
			exit()
		ret = process_line(line, consequence_header_fields)
		count += 1
		if count%100000==0: print "%dK lines out of 15014K (%5.2f%%)" % (count/1000,  float(count)/15014744*100)
		if not ret: continue
		[chrom, addr, ref, variants, consequences,  ac, an, ac_afr, an_afr,
		 ac_amr, an_amr, ac_asj, an_asj, ac_eas, an_eas, ac_fin, an_fin,
		 ac_nfe, an_nfe, ac_oth, an_oth, ac_sas, an_sas] = ret
		if False:
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
		#store_or_update(cursor, table, fixed_fields, update_fields, verbose=False, primary_key='position')
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914