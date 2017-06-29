#!/usr/bin/python

from integrator_utils.mysql import *

#########################################
def parse_csq(csq_string):
	fields = csq_string.split(",")
	genes_and_names = set([  tuple(field.split("|")[3:5])   for field in fields])
	print genes_and_names
	exit()
	return

#########################################
def parse_info(info):
	fields = info.split(";")
	named_fields= dict (map(lambda x: x.split("="), [f for f in fields if '=' in f]))
	#genes_affected = parse_csq(named_fields['CSQ']) #, actually I can live without that info here
	return named_fields

#########################################
def process_line(line):
	fields = line.split("\t")
	[chrom, addr, junk, ref, variants, quality_score,  filter, info_string] = fields[:8]
	if filter!='PASS': return None
	info = parse_info (info_string)
	if info["AC_Adj"]=='0' : return None
	return_list  = [chrom, addr, ref, variants]
	return_list += [info["AC_Adj"],  info["AN_Adj"]]
	return_list += [info["AC_AFR"],  info["AN_AFR"]]
	return_list += [info["AC_AMR"],  info["AN_AMR"]]
	return_list += [info["AC_EAS"],  info["AN_EAS"]]
	return_list += [info["AC_FIN"],  info["AN_FIN"]]
	return_list += [info["AC_NFE"],  info["AN_NFE"]]
	return_list += [info["AC_OTH"],  info["AN_OTH"]]
	return_list += [info["AC_SAS"],  info["AN_SAS"]]

	return return_list


##########################################
def main():
	infile = open("/databases/exac/ExAC_nonTCGA.r1.sites.vep.vcf")

	db, cursor = connect()

	reading = False
	for line in infile:
		if not reading:
			if line[:6] == '#CHROM': reading=True
			continue
		ret = process_line(line)
		if not ret: continue
		[chrom, addr, ref, variants, ac, an, ac_afr, an_afr,
		 ac_amr, an_amr, ac_eas, an_eas, ac_fin, an_fin,
		 ac_nfe, an_nfe, ac_oth, an_oth, ac_sas, an_sas] = ret
		table = "exac_freqs_chr_" + chrom
		fixed_fields  = {'position':int(addr)}
		update_fields = {'reference':ref, 'variants':variants, 'variant_counts':ac, 'total_count':an}
		update_fields['afr_counts'] = ac_afr
		update_fields['afr_tot_count'] = an_afr
		update_fields['amr_counts'] = ac_amr
		update_fields['amr_tot_count'] = an_amr
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