#!/usr/bin/python3

from integrator_utils.python.mysql import *
from integrator_utils.python.processes import parallelize
from random import shuffle

# output the tables rather than trying to store them from python - much faster
# but do the sanity checking in the following script first
'''
i=1
sudo mysqlimport  --local gnomad freqs_chr_$i.tab
to get prompted for password
then

 for i in $(seq 2 22)
> do
> echo $i
> sudo mysqlimport  --local gnomad freqs_chr_$i.tab
> done
repeat for X and Y
for i in X Y; do echo $i; done

"NOTE that the file name must match the table name:
For each text file named on the command line, mysqlimport strips any extension 
from the file name and uses the result to determine the name of the table into 
which to import the file's contents. For example, files named patient.txt,
patient.text , and patient all would be imported into a table named patient.
Also before you begin, do
mysql> SET GLOBAL local_infile = 1;
"

'''
# used to be ["AC", "AN", "AF", "GC", "Hemi"]
count_fields = ["AC", "AN", "AF"]
# used to be uppercase, now is lower
population_fields = [field.lower() for field in ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE","OTH","SAS"]]
# now we also have subpopulations  as in
# AF_nfe_bgr        Alternate allele frequency in samples of Bulgarian ancestry
# AF_nfe_est        Alternate allele frequency in samples of Estonian ancestry
# but we are going to ignore that here

gender = [field.lower() for field in ["Male","Female"]]

all_keys = ['chrom', 'addr', 'ref', 'variants',   'ac', 'an', 'ac_afr', 'an_afr',
            'ac_amr', 'an_amr', 'ac_asj', 'an_asj', 'ac_eas', 'an_eas', 'ac_fin', 'an_fin',
			'ac_nfe', 'an_nfe', 'ac_oth', 'an_oth', 'ac_sas', 'an_sas', 'nhomalt']
counts = [key for key in all_keys if key[:2]=='ac']
totals = [key for key in all_keys if key[:2]=='an']
specials = ['nhomalt']
# the annotation is crap ...
# 22:17669352 G>GGGACA is an SNV in somebody's book

verbose = False


##########################################
def format_tsv_line(named_fields):

	# not sure why this happens - filtered varaints (the ones that do not have 'PASS')
	# still may have count > 0
	# for the gnomad blogpost
	# "All sites with no high quality genotype (AC = 0) are marked as filtered using the AC0 filter."
	if named_fields['ac']== '0':
		print(named_fields)
		exit()
		return

	# in trying for a compact storage, gnomad packs several variants into a single line
	# for example ref CT, variants C,AT to account for a deletion and an SNV
	# we are going to turn that back int CT, C  and C, A

	fixed_fields  = {'position':int(named_fields['addr']),  'reference': named_fields['ref'], 'variant':named_fields['variants']}
	update_fields = {}

	for count in counts:
		if count == 'ac':
			new_key = 'variant_count'
		else:
			population = count.split("_")[-1]
			new_key = "_".join([population,'count'])
		update_fields[new_key] = int(named_fields[count])

	for count in totals:
		if count == 'an':
			new_key = 'total_count'
		else:
			population = count.split("_")[-1]
			new_key = "_".join([population,'tot','count'])
		update_fields[new_key] = int(named_fields[count])

	for count in specials:
		if count == 'nhomalt':
			new_key = 'homozygote_count'
		else:
			new_key = count
		update_fields[new_key] = int(named_fields[count])

	all_fields = fixed_fields
	all_fields.update(update_fields)
	# getting things in the right order
	db_field_names = ['position', 'reference', 'variant', 'variant_count', 'total_count',
						'afr_count', 'afr_tot_count', 'amr_count', 'amr_tot_count', 'asj_count',
						'asj_tot_count', 'eas_count', 'eas_tot_count', 'fin_count', 'fin_tot_count',
						'nfe_count', 'nfe_tot_count', 'oth_count', 'oth_tot_count', 'sas_count',
						'sas_tot_count', 'homozygote_count', 'hotspot_id']
	all_fields['hotspot_id'] = "\\N"
	tabbed_line = "\t".join([str(all_fields[field_name]) for field_name in db_field_names])
	return tabbed_line

#########################################
def parse_counts(chrom,  info):

	# AC        Alternate allele count for samples
	# AF        Alternate allele frequency in samples
	# AN        Total number of alleles in called genotypes
	fields = info.split(";")
	named_fields = dict ([x.split("=") for x in [f for f in fields if '=' in f]])
	data = {}
	for cnt in count_fields:
		data[cnt] = {}
		for p in population_fields+['overall']:
			data[cnt][p] = {}
			for g in gender+['both']: data[cnt][p][g] = '0'

	# parsing field names of the type 'AF_eas_male' (for example)
	homozygotes = 0
	for k, v in named_fields.items():
		if k=="nhomalt":
			homozygotes = v
			continue
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
			elif field[1] in ['raw','popmax']:
				continue
			else:
				print("unrecognized key:", k)
				exit()

		elif l==3:
			[p,g] = field[1:3]
			if not p in population_fields or not g in gender:
				if verbose: print("unrecognized key:", k, "(subpopulation?)")
				continue
			data[count_type][p][g] = v

		else:
			print("unrecognized key:", k)
			exit()



	# sanity check
	if verbose:
		if chrom in ['X','Y']:
			for p in population_fields:

				for g in gender+['both']:
					print(p, g, data['AC'][p][g], data['AN'][p][g], data['AF'][p][g])
			print('overall', data['AC']['overall']['both'], data['AN']['overall']['both'], data['AF']['overall']['both'])
			print('-----')
		else:
			# AN = total number of alleles
			total = 0
			for p in population_fields:
				total += int(data['AN'][p]['both'])
				print(p, data['AC'][p]['both'], data['AN'][p]['both'], data['AF'][p]['both'])
			print(data['AC']['overall']['both'], data['AN']['overall']['both'], data['AF']['overall']['both'])
			if int(data['AN']['overall']['both']) != total:
				print("count mismatch:",data['AN']['overall']['both'], total)
				exit()

	counts = []
	counts.append(data['AC']['overall']['both'])
	counts.append(data['AN']['overall']['both'])
	for p in population_fields:
		counts.append(data['AC'][p]['both'])
		counts.append(data['AN'][p]['both'])
	counts.append(homozygotes)

	return counts


########################################
def parse_data_line(line):

	fields = line.rstrip().split("\t")
	[chrom, addr, junk, ref, variant, quality_score,  filter, info_string] = fields[:8]
	if filter!='PASS': return None
	counts = parse_counts (chrom, info_string)

	return [chrom, addr, ref, variant] + counts


##########################################
def process_line(line,):

	ret = parse_data_line(line)

	if not ret: return  # no PASS line

	# unpack variants
	named_fields = dict(list(zip(all_keys, ret)))
	# in trying for a compact storage, gnomad 2.0.1 packed several variants into a single line
	# for example ref CT, variants C,AT to account for a deletion and an SNV
	# mercifully, it seems we do not have that piece of idiocy present any more
	return format_tsv_line(named_fields)


##########################################
def process_chromosome_data(chromosomes, other):
	[indir, outdir, mysql_conf_file] = other

	db     = connect_to_mysql(conf_file=mysql_conf_file)
	cursor = db.cursor()

	for chrom in chromosomes:
		infnm  = "%s/gnomad.exomes.r2.1.1.sites.%s.liftover_grch38.vcf" % (indir, chrom)
		outfnm = "%s/freqs_chr_%s.tab" % (outdir, chrom)

		infile  = open(infnm,"r")
		outfile = open(outfnm,"w")

		reading = False
		count = 0

		for line in infile:
			if not reading:
				if line[:6] == '#CHROM': reading=True
				continue
			chrom = line.split("\t")[0]
			count += 1
			if count%10000==0: print("%dK lines of  %s" % (count/1000,  chrom))
			tabbed_line = process_line(line)
			if not tabbed_line: continue # no PASS
			outfile.write("%d\t%s\n"%(count, tabbed_line))

		infile.close()
		outfile.close()

	cursor.close()
	db.close()

##########################################
def main():

	# This one is good. From gnomad blogpost, literally:
	# https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad/
	# "There is no chromosome Y coverage / calls in the gnomAD genomes, because reasons."
	# this seems to have been fixed in gnomad 2.1.1
	#chromosomes = [str(i) for i in range(1,23)] + ['X', 'Y']
	indir = "/storage/databases/gnomad"
	outdir = "/home/ivana/scratch/gnomad_freqs_tmp"
	mysql_conf_file = "/home/ivana/.tcga_conf"
	for dep in [indir, outdir, mysql_conf_file]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	number_of_chunks = 1
	chromosomes = ['12']
	# number_of_chunks = 8
	# chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
	# shuffle(chromosomes) # shuffles in place and returns None
	parallelize(number_of_chunks, process_chromosome_data, chromosomes, [indir, outdir, mysql_conf_file])


	return


#########################################
if __name__ == '__main__':
	main()

# 8,608,914