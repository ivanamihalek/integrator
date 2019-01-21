#!/usr/bin/python


from integrator_utils.python.mysql import *
from random import randint
from numpy.random import normal
def hdr(member):
	return "##fileformat=VCFv4.1\n" \
			"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" \
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n" \
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" \
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT %s\n" % member


##########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"gnomad")
	member = "father"
	outf = open ("%s.vcf"%member,"w")
	outf.write(hdr(member))
	tot = 0
	selected = 0
	for chrom in [str(i) for i in range(1,23)] + ['X']: # no 'Y' in gnomadwc -l
		table = "gnomad_freqs_chr_" + chrom
		print table
		qry = "select position, reference, variant, variant_count, total_count from %s" % table
		ret = search_db(cursor, qry)
		if not ret:
			print "no ret for", qry
			exit(1)
		for [position, reference, variant, variant_count, total_count] in ret:
			tot += 1
			# random.randint(a, b):  Return a random integer N such that a <= N <= b.
			if randint(1,total_count)>variant_count: continue
			# I'm not sure what should be the best thing to have here, so for simulation purposes, I will take  the
			# the number of quality bases (DP) to be equal to the sum of allelic depths
			dp  = int(normal(100,50,1)) # norma(center, width, sample size)
			if dp<4: continue
			ad1 = int(normal(dp/2,dp/10,1))
			if ad1<0: ad1=0
			if ad1>dp: ad1=dp
			if ad1<dp/4:
				gt="0/1"
			elif ad1>3*dp/4:
				gt="1/0"
				continue # this is a non-variant
			else:
				gt="1/1"

			outf.write("\t".join(["chr%s"%chrom, str(position), ".", reference, variant,
						str(randint(100,5000)), ".",  ".", "GT:AD:DP",
						"{}:{},{}:{}".format(gt, ad1, dp-ad1, dp)])+"\n")
			selected += 1

	print "selected {} out of {}".format(selected, tot)
	outf.close()
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
