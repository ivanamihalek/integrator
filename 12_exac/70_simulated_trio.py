#!/usr/bin/python3

# gnomad is based on hg19 (Jan 2019)
# http://gnomad.broadinstitute.org/faq

from integrator_utils.python.mysql import *
from integrator_utils.python.restfuls import ucsc_fragment_sequence
from random import randint, choice
from numpy.random import normal
def hdr(member):
	return "##fileformat=VCFv4.1\n" \
			"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" \
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n" \
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" \
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT %s\n" % member

# https://app.face2gene.com, chrome
##########################################
def cookup_depths(skew="none"):
	dp = 0
	ad1 = 0
	gt="1/0"

	# I'm not sure what should be the best thing to have here,
	# so for simulation purposes, I will take  the
	# the number of quality bases (DP) to be equal to the sum of allelic depths
	dp  = int(normal(100,50,1)) # norma(center, width, sample size)
	if dp<4: return[dp, ad1, gt]

	ad1 = int(normal(dp/2,dp/10,1))
	if skew=="left":
		ad1 += dp/2
	elif skew=="right":
		ad1 -= dp/2

	if ad1<0: ad1=0
	if ad1>dp: ad1=dp

	if ad1<dp/4:
		gt="1/1"  # variant on both alleles
	elif ad1>3*dp/4:
		gt="0/0"  # reference on both alleles
	else:
		gt="0/1"  # heterozygous
	return [dp,ad1,gt]


##########################################
def gnomad_variants(cursor, chrom, outf):
	table = "gnomad_freqs_chr_" + chrom
	print(table)
	qry = "select position, reference, variant, variant_count, total_count from %s" % table
	# IT LOOKS LIKE I HAVA BUG IN HTE WAY I AM FILLING TH GNOMAD, SEE chr 22 16279328
	# ref and variant are th same , while this should actually be a deletion
	ret = search_db(cursor, qry)
	if not ret:
		print("no ret for", qry)
		exit(1)
	for [position, reference, variant, variant_count, total_count] in ret:
		if reference==variant: continue # see the bug described above
		alleles = {}
		for parent in ['mother', 'father']:
			alleles[parent] = [reference,reference]
			# random.randint(a, b):  Return a random integer N such that a <= N <= b.
			if randint(1,total_count)>variant_count: continue # we didn't pick up a rare variant
			[dp, ad1, gt] = cookup_depths()
			if dp<4: continue # we didn't pass the quality check
			if gt=="0/0": continue # this is a non-variant
			if gt=="1/1":
				alleles[parent] = [variant, variant]
			else:
				alleles[parent] = [reference, variant]
			outf[parent].write("\t".join(["chr%s"%chrom, str(position), ".", reference, variant,
								str(randint(100,5000)), ".",  ".", "GT:AD:DP",
								 "{}:{},{}:{}".format(gt, ad1, dp-ad1, dp)])+"\n")

		# did the child inherit the variant (I get in trouble here if someone actually knows the haplotype)
		if variant in alleles["mother"] or variant in alleles["father"]:
			alleles["child"] = [ alleles["mother"][randint(0,1)], alleles["father"][randint(0,1)] ]
			if variant in alleles["child"]:
				skew = "right" if alleles["child"][0] == alleles["child"][1] else "none"
				[dp, ad1, gt] = cookup_depths(skew)
				outf["child"].write("\t".join(["chr%s"%chrom, str(position), ".", reference, variant,
								str(randint(100,5000)), ".",  ".", "GT:AD:DP",
								"{}:{},{}:{}".format(gt, ad1, dp-ad1, dp)])+"\n")

##########################################
def random_variants(cursor, assembly, chrom, member, outf):
	nucleotides = {'A','C','T','G'}
	table = "ucsc.refgenes"
	qry = "select distinct(name) from %s where chromosome='chr%s'" % (table,chrom)
	names = [line[0] for line in  search_db(cursor, qry)]
	for name in names:
		if (randint(1,20)>2): continue
		qry = "select exon_count, exon_starts, exon_ends from %s where chromosome='chr%s' and name='%s' limit 1" % (table, chrom, name)
		[exon_count, exon_sts, exon_es] = search_db(cursor, qry)[0]
		exon_starts = [int(i) for i in exon_sts.decode().split(",")] # withtout decode() get  a bytes-like object is required, not 'str'
		exon_ends   = [int(i) for i in exon_es.decode().split(",")]
		random_exon = randint(1,exon_count)-1
		random_pos  = randint(exon_starts[random_exon], exon_ends[random_exon])
		nucleotide  = ucsc_fragment_sequence(assembly, chrom, random_pos, random_pos).upper()
		if (randint(1,10)<2):
			mutation = nucleotide+nucleotide
		else:
			mutation = choice(list(nucleotides-{nucleotide}))
		[dp, ad1, gt] = cookup_depths()
		if dp<4: continue # we didn't pass the quality check
		if gt=="0/0": continue # this is a non-variant
		outf[member].write("\t".join(["chr%s"%chrom, str(random_pos), ".", nucleotide, mutation,
								str(randint(100,5000)), ".",  ".", "GT:AD:DP",
								"{}:{},{}:{}".format(gt, ad1, dp-ad1, dp)])+"\n")
		if member in ["mother","father"] and randint(1,4)<2:
			[dp, ad1, gt] = cookup_depths()
			if dp<4: continue # we didn't pass the quality check
			if gt=="0/0": continue # this is a non-variant
			outf["child"].write("\t".join(["chr%s"%chrom, str(random_pos), ".", nucleotide, mutation,
								str(randint(100,5000)), ".",  ".", "GT:AD:DP",
								"{}:{},{}:{}".format(gt, ad1, dp-ad1, dp)])+"\n")


##########################################
def main():

	assembly = 'hg19'

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"gnomad")
	outf = {}
	for member in ["mother", "father","child"]:
		outf[member] = open ("%s.vcf"%member,"w")
		outf[member].write(hdr(member))
	tot = 0
	selected = 0
	# at the moment I am not interested in sex-chromosome mutations
	#for chrom in [str(i) for i in range(1,23)]: #+ ['X']: # no 'Y' in gnomad
	for chrom in ['22']:
		# gnomad variants
		gnomad_variants(cursor, chrom, outf)
		# idiosyncratic/de novo variants in parents
		random_variants(cursor, assembly,  chrom, "mother", outf)
		random_variants(cursor, assembly,  chrom, "father", outf)
		# idiosyncratic/de novo variants in child
		random_variants(cursor, assembly,  chrom, "child", outf)

	# disease variant(s) in the child

	for member in ["mother", "father","child"]:
		outf[member].close()
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
