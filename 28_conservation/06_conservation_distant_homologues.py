#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess
from math import log
from Bio.Align.Applications import MuscleCommandline

scratch      = "/home/ivana/scratch"
blastp       = "/usr/local/bin/blastp"
#uniprotdb    = "/databases/uniprot/blast/uniprot_sprot.fasta"
# run on bronto to use trembl
uniprotdb    = "/databases/uniprot/blast/uniprot_trembl.fasta"
blastextract = "/usr/local/bin/blastdbcmd"
afa2msf      = "/home/ivana/pypeworks/integrator/integrator_utils/afa2msf.pl"
restrict     = "/home/ivana/pypeworks/integrator/integrator_utils/restrict_msf_to_query.pl"
msf2afa      = "/home/ivana/pypeworks/integrator/integrator_utils/msf2afa.pl"
structure_repo = "/home/ivana/monogenic/public/pdb"

if gethostname()=='brontosaurus':
	mono_db   = "monogenic_production"
	blimps_db = "blimps_production"
else:
	mono_db   = "monogenic_development"
	blimps_db = "blimps_development"

##########################################
def column_entropy(string):
	bin = {}
	for char in string:
		if not bin.has_key(char): bin[char]=0
		bin[char] +=1
	entropy = 0.0
	for char, count in bin.iteritems():
		freq = float(count)/len(string)
		if count>0:
			entropy -= freq*log(freq)
	return entropy/log(20)

##########################################
def blastsearch(sequence,uniprot_id):

	queryfile = "{}/{}.qry.fasta".format(scratch,uniprot_id)
	outfile = open(queryfile,"w")
	outfile.write(">{}_query\n{}\n".format(uniprot_id, sequence))
	outfile.close()
	#blast search against uniprot
	outfile =  "{}/{}.blastout".format(scratch,uniprot_id)
	if not os.path.exists(outfile) or os.stat(outfile).st_size == 0:
		cmd_format = "{} -db {} -query {} -out {} -evalue 1.0e-20  -outfmt 6 -max_target_seqs 1000"
		cmd = cmd_format.format(blastp, uniprotdb, queryfile, outfile)
		subprocess.call(cmd, shell=True)
	lowest_e = subprocess.check_output("tail -n 1 {}".format(outfile), shell=True).split()[-2]
	# find the corresponding seqs
	fastafile = "{}/{}.fasta".format(scratch,uniprot_id)
	if not os.path.exists(fastafile) or os.stat(fastafile).st_size == 0:
		idfile    = "{}/{}.ids".format(scratch,uniprot_id)
		cmd = "awk '{print $2}' %s > %s" %(outfile,idfile)
		subprocess.call(cmd, shell=True)
		cmd = "{} -db {} -entry_batch {} -out {}".format(blastextract, uniprotdb, idfile, fastafile)
		subprocess.call(cmd, shell=True)
		# prepend query sequence - moreutils must be installed
		cmd = "cat {} {} | sponge {}".format(queryfile, fastafile, fastafile)
		subprocess.call(cmd, shell=True)
	return fastafile, lowest_e

##########################################
def align(fastafile, uniprot_id):
	afafile = "{}/{}.afa".format(scratch,uniprot_id)
	if not os.path.exists(afafile) or os.stat(afafile).st_size == 0:
		cline = MuscleCommandline(input=fastafile, out=afafile,  verbose=True)
		stdout, stderr = cline()
	return afafile

##########################################
def conservation(restricted_msf, uniprot_id):
	# find structure file
	outf = open("{}.cmd".format(uniprot_id),"w")
	outf.write("sink 0.3\nskip_qry\nmethod rvet\n")
	outf.write("align {}\nrefseq {}_query\n".format(restricted_msf,uniprot_id))
	outf.write("outn {}\n".format(uniprot_id))
	outf.close()
	exit()
	return None

##########################################
def conserved_columns(restricted_afa, uniprot_id):
	infile = open(restricted_afa,"r")
	seq = {}
	for line in  infile:
		line = line.rstrip()
		if len(line)==0: continue
		if line[0]=='>':
			seqname = line[1:]
		else:
			if not seq.has_key(seqname): seq[seqname] =""
			seq[seqname] += line
	infile.close()

	columns = ['']*len(seq[uniprot_id+"_query"])
	for name, sequence in seq.iteritems():
		for i in range(len(sequence)): columns[i]+=sequence[i]
	conservation = []
	for column in columns:
		conservation.append(1-column_entropy(column))
	return conservation

##########################################
def main():

	for dependency in [scratch, blastp, uniprotdb, blastextract, afa2msf, restrict, msf2afa]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()
	db, cursor = connect()
	qry = "select * from %s.diseases" % mono_db
	ret = search_db(cursor, qry)
	plotdata = ""
	for line in ret:
		print qry
		print ret
		exit()
		[id, name_short, name_long, omim_ids, description, prim, sec] = line
		print "="*60
		print name_short
		for omim_id in omim_ids.split(";"):
			qry = "select approved_symbol, ensembl_gene_id from omim_genemaps where mim_number='%s'" % omim_id
			ret2 = search_db(cursor, qry)
			for line2 in ret2:
				approved_symbol = line2[0]
				# ensembl_gene_id = line2[1]
				qry = "select uniprot_id from uniprot_basic_infos where gene_name='%s'" % approved_symbol
				ret3 = search_db(cursor, qry)
				if not ret3:
					print "uniprot id not found"
					exit()
				uniprot_id = ret3[0][0]
				print uniprot_id
				qry = "select sequence from %s.uniprot_seqs where uniprot_id='%s'" % (mono_db,uniprot_id)
				ret3 = search_db(cursor, qry)
				if not ret3:
					print "sequence  not found"
					exit()
				sequence = ret3[0][0]
				#qry = "select structure_file from structures where gene_symbol='%s'" %  approved_symbol
				#ret3 = search_db(cursor, qry)
				#if not ret3:
				#	print "structure  not found"
				#	continue
				#structure_file = ret3[0][0]
				#blasr
				fastafile, lowest_e  = blastsearch(sequence,uniprot_id)
				afafile   = align(fastafile,uniprot_id)
				# afa2ms
				msffile = afafile.replace("afa","msf")
				subprocess.call("{} {} > {}".format(afa2msf,afafile,msffile), shell=True)
				# restrict
				restricted_msf = msffile.replace("msf","restr.msf")
				subprocess.call("{} {} {}_query > {}".format(restrict,msffile, uniprot_id, restricted_msf), shell=True)
				# conservation
				#score_file = conservation(restricted_msf, uniprot_id)
				# cluster?
				restricted_afa = afafile.replace("afa","restr.afa")
				subprocess.call("{} {} > {}".format(msf2afa,restricted_msf,restricted_afa), shell=True)

				position_conservation = conserved_columns(restricted_afa, uniprot_id)
				conserved_positions = []
				for i in range(len(position_conservation)):
					if position_conservation[i] ==1:  conserved_positions.append(str(i+1))

				cons_string = ",".join(conserved_positions)
				print cons_string
				plotdata += "%5.1f  %5.1f   %s  %d \n"%\
				            (-log(float(lowest_e),10),len(cons_string)*100.0/len(sequence),name_short,len(sequence))
				#qry  = "update monogenic_development.uniprot_seqs "
				#qry += "set conservation_in_vertebrates='%s' " % cons_string
				#qry += "where uniprot_id='%s' " % uniprot_id
				#search_db(cursor, qry, verbose=False)

	print plotdata

########################################
if __name__ == '__main__':
	main()
