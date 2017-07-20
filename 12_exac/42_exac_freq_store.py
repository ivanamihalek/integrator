#!/usr/bin/python

from integrator_utils.mysql import *

count_fields = ["AC", "AN", "AF", "GC", "Hemi"]
population_fields = ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE","OTH","SAS"]
gender = ["Male","Female"]

#########################################
def parse_csq(csq_string):
	fields = csq_string.split(",")
	genes_and_names = set([tuple(field.split("|")[3:5])   for field in fields])
	print genes_and_names
	exit()
	return

#########################################
def parse_info(chrom, info):

	fields = info.split(";")
	named_fields= dict (map(lambda x: x.split("="), [f for f in fields if '=' in f]))

	data = {}
	for cnt in count_fields:
		data[cnt] = {}
		for p in population_fields+['overall']:
			data[cnt][p] = {}
			for g in gender+['both']: data[cnt][p][g] = 0

	for  k,v  in named_fields.iteritems():
		field = k.split("_")
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
	return counts

#########################################
def process_line(line):
	fields = line.split("\t")
	[chrom, addr, junk, ref, variants, quality_score,  filter, info_string] = fields[:8]
	if filter!='PASS': return None
	counts = parse_info (chrom, info_string)


	return [chrom, addr, ref, variants] + counts


##########################################
def main():
	infile = open("/databases/exac/gnomad.exomes.r2.0.1.sites.vcf")
	#infile = open("/databases/exac/testY.vcf")

	db, cursor = connect()

	reading = False
	count = 0
	for line in infile:
		if not reading:
			if line[:6] == '#CHROM': reading=True
			continue
		ret = process_line(line)
		count += 1
		if count%100000==0: print "%dK lines out of 15014K (%5.2f%%)" % (count/1000,  float(count)/15014744*100)
		if not ret: continue
		[chrom, addr, ref, variants, ac, an, ac_afr, an_afr,
		 ac_amr, an_amr, ac_asj, an_asj, ac_eas, an_eas, ac_fin, an_fin,
		 ac_nfe, an_nfe, ac_oth, an_oth, ac_sas, an_sas] = ret
		table = "exac_freqs_chr_" + chrom
		fixed_fields  = {'position':int(addr)}
		update_fields = {'reference':ref, 'variants':variants, 'variant_counts':ac, 'total_count':an}
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