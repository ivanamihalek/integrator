#!/usr/bin/python

import shlex

##########################################
def main():
	infile = open("/databases/exac/ExAC_nonTCGA.r1.sites.vep.vcf")
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
		print fields[0].split('=')[-1], "      ", fields[-1].split('=')[-1]

	print
	print '''It looks like the most interesting piece of info are the adjusted allele counts:

=========================================================
Population independent annotation (also included in VCF header):

Adjusted Alt Allele Counts (DP >= 10 & GQ >= 20)
##INFO=<ID=AC_Adj,Number=A,Type=Integer,Description="Adjusted Allele Counts">
AC_Adj <= AC

Number of Heterozygous Individuals (DP >= 10 & GQ >= 20)
##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Adjusted Heterozygous Counts">

Number of Homozygous Alt Allele Individuals (DP >= 10 & GQ >= 20)
##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Adjusted Homozygous Counts">

For chr1-22:
sum(AC_Adj) = sum(AC_Het) + 2*sum(AC_Hom)

=========================================================
Adjustments made on chrX and chrY

Number of Hemizygous Alt Allele Individuals (DP >= 10 & GQ >= 20) Note: ONLY appears on chrX (non-PAR) and chrY
##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Adjusted Hemizygous Counts">

AC_Hemi is a count Male alt alleles, where each site (excluding PAR) only has one allele.
AC_Hom is a count of Female alt alleles
AC_Het is a count of Female alt alleles

For chrX (non-PAR)
sum(AC_Adj) = sum(AC_Hemi) + sum(AC_Het) + 2*sum(AC_Hom)
AN_Adj (and all population AN on chrX) = 2*n_Female + n_Male

For chrY
sum(AC_Adj) = sum(AC_Hemi)
AN_Adj (and all population AN on chrY) = n_Male

Pseudoautosomal regions (PARs) were taken from http://genome.ucsc.edu/cgi-bin/hgGateway
X:60001-2699520
X:154931044-155260560

=========================================================

By comparison with the website it seems that the number that they cite as the allele frequency is AC_Adj/AN_Adj
there is AF or allele frq number in the VCF, but not clear what that is.

Whrn comparing with the ExAC browser keep in mind that the browser itself might not be completely up to date.

	'''

	return
#########################################
if __name__ == '__main__':
	main()
