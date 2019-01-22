#!/usr/bin/python3

# gnomad is based on hg19 (Jan 2019)
# http://gnomad.broadinstitute.org/faq

from integrator_utils.python.mysql import *
from integrator_utils.python.restfuls import ucsc_fragment_sequence
from random import randint, choice
from numpy.random import normal
import subprocess

##########################################
def read_variants(infile):
	vars = {}
	zygosity = {}
	for  member in ["mother", "father","child"]:
		zygosity[member] = {}
		inf = open(infile[member], "r")
		for line in inf:
			fields = line.split()
			if not "chr" in fields[0][0:3]: continue
			chrom = fields[0][3:]
			position_from = int(fields[1])
			ref = fields[3]
			var = fields[4]
			position_to = position_from + len(ref) -1
			if not chrom in vars: vars[chrom] = []
			variant = " ".join([chrom, str(position_from), str(position_to), ref, var])
			if not variant in vars[chrom]: vars[chrom].append(variant)
			zygosity[member][variant] = fields[-1].split(":")[0]
		inf.close()
	return vars, zygosity

##########################################
def annovar_input(vars, outfile):
	outf = open(outfile, "w")
	for chrom,vars in vars.items():
		for var in vars:  outf.write(var+"\n")
	outf.close()

##########################################
def annovar_annotation(assembly, vars):
	base_name = "test"
	avinfile = "test.avinput"
	annovar_input(vars,avinfile)

	avoutname = "%s.%s_multianno.txt" % (base_name, assembly)
	cmd  = "/home/ivana/third/annovar/table_annovar.pl %s " % avinfile
	cmd += "/home/ivana/third/annovar/humandb/ -buildver %s -out %s " % (assembly, base_name)
	cmd += " -protocol refGene  -operation g  -nastring . "
	print(cmd)
	subprocess.call(cmd, shell=True)
	# clean the junk
	cmd = "rm %s.refGene.variant_function " % base_name
	cmd +="%s.refGene.exonic_variant_function %s.refGene.log" % (base_name, base_name)
	subprocess.call ( cmd, shell=True)

	annotation = {}
	inf = open(avoutname, "r")
	for line in inf:
		fields = line.rstrip().split("\t")
		if fields[0]=="Chr": continue
		variant = " ".join(fields[:5])
		annotation[variant] = fields[5:7] + fields[8:]

	return annotation

##########################################
def main():

	assembly = 'hg19'
	disease = "ANGELMAN"
	db = connect_to_mysql()
	infile = {}
	for member in ["mother", "father","child"]:
		infile[member] = "%s_%s.vcf"%(disease[0].lower(),member)

	outfile = "%s_%s.vcf"%(disease[0].lower(),"trio")
	outf = open(outfile, "w")
	outf.write("\t".join(["chrom", "position", "ref", "alt", "zygosity mom", "zygosity dad",  "zygosity proband",
						"location", "gene name", "effect","protein change"]) + "\n")
	vars, zygosity = read_variants(infile)
	annotation = annovar_annotation(assembly, vars)
	for var,annotation_fields in annotation.items():
		[chrom, position_from, position_to, ref, alt] = var.split(" ")
		zyg_descriptive = []
		for member in ["mother", "father", "child"]:
			if var in zygosity[member]:
				z = zygosity[member][var]
				if z=="0/0":
					zyg_descriptive.append("no_variant")
				elif z=="0/1":
					zyg_descriptive.append("heterozygous")
				elif z=="1/1":
					zyg_descriptive.append("homozygous")
			else:
				zyg_descriptive.append("no_variant")
		outf.write("\t".join([chrom, position_from, ref, alt] + zyg_descriptive + annotation_fields) +"\n")

	outf.close()

	return True

#########################################
if __name__ == '__main__':
	main()
