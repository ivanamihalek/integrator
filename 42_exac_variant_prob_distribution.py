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
	info = parse_info (info_string)
	if info["AC_Adj"]=='0' : return None
	return  [chrom, addr, ref, variants, info["AC_Adj"],  info["AN_Adj"]]


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
		[chrom, addr, ref, variants, ac,  an] = ret
		table = "exac_freqs_chr_" + chrom
		fixed_fields  = {'position':int(addr)}
		update_fields = {'reference':ref, 'variants':variants, 'variant_counts':ac, 'total_count':an}
		store_or_update(cursor, table, fixed_fields, update_fields, verbose=True, primary_key='position')
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914