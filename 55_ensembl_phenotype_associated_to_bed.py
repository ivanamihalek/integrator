#!/usr/bin/python


import re

##########################################
# Homo_sapiens_phenotype_associated.vc columns:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
def main():

	infile = open("/databases/ensembl/hg19/Homo_sapiens_phenotype_associated.vcf","r")

	for line in infile:
		if line[0]=='#': continue
		fields =  line.rstrip().split("\t")
		if not len(fields)==8: continue
		if "benign" in fields[-1].lower(): continue
		descr = ";".join([fields[3], fields[4], fields[-1]])
		pos = int(fields[1])
		start = pos-10
		stop  = pos+10
		chrom = fields[0]
		print "\t".join([chrom, str(start), str(stop), descr])
	return True

#########################################
if __name__ == '__main__':
	main()
