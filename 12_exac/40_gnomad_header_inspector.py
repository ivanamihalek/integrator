#!/usr/bin/python

import shlex

##########################################
def main():
	#infile = open("ExAC_nonTCGA.r1.sites.vep.vcf")
	infile = open("/databases/exac/gnomad_headers.txt")
	for line in infile:
		if line[:len("##reference")] == "##reference":
			print "reference assembly:",line.rstrip().split("=")[-1]
			continue
		if line[:6] == '#CHROM': break
		if line[:len('##INFO=<')] != '##INFO=<': continue
		my_splitter = shlex.shlex(line.rstrip()[len('##INFO=<'):-1], posix=True)
		my_splitter.whitespace += ','
		my_splitter.whitespace_split = True
		fields = list(my_splitter)
		field_name =  fields[0].split('=')[-1]
		field_content = fields[-1].split('=')[-1]
		if field_name=='CSQ':
			print field_name,field_content.split(':')[0]
			for subcontent in field_content.split(':')[1].split('|'):
				print "\t", subcontent
		else:
			print  field_name, "      ",  field_content

	print

	return
#########################################
if __name__ == '__main__':
	main()
