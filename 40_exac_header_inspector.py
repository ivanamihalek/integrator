#!/usr/bin/python

import shlex

##########################################
def main():
	infile = open("/databases/exac/ExAC_nonTCGA.r1.sites.vep.vcf")
	for line in infile:
		if line[:6] == '#CHROM': break
		if line[:len('##INFO=<')] != '##INFO=<': continue
		my_splitter = shlex.shlex(line.rstrip()[len('##INFO=<'):-1], posix=True)
		my_splitter.whitespace += ','
		my_splitter.whitespace_split = True
		fields = list(my_splitter)
		print fields[0].split('=')[-1], "      ", fields[-1].split('=')[-1]

	return
#########################################
if __name__ == '__main__':
	main()
