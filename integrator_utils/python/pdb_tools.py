
import pdb_constants as pdbc
import os, subprocess, shutil
from collections import OrderedDict

blastp       = "/usr/local/bin/blastp"
pdb_blast_db = "/databases/pdb/blast/pdb_seqres.fasta" # expected to be formatted using makeblastdb
pdb_name_resolution_table = "/databases/pdb/blast/name_resolution.txt"


##########################################
def pdb_check_resources():
	for rsrc in [blastp, pdb_blast_db, pdb_name_resolution_table]:
		if not os.path.exists(rsrc):
			print rsrc, "not found"
			return False
		if  os.path.getsize(rsrc)==0:
			print rsrc, "is empty"
			return False
	return True


##########################################
def pdb_to_sequence(path, filename):
	sequence = {}
	infile = open (path+"/"+filename, "r")
	for line in infile:
		if line[:4] != 'ATOM': continue
		chain = line[pdbc.chain_pos]
		resname = line[pdbc.res_name_pos:pdbc.res_name_pos+pdbc.res_name_length]
		resnumber = int(line[pdbc.res_number_pos:pdbc.res_number_pos+pdbc.res_number_length])
		if not pdbc.aa_translation.has_key(resname): continue
		if not sequence.has_key(chain): sequence[chain] = {}
		sequence[chain][resnumber] = pdbc.aa_translation[resname]

	return sequence


##########################################
def resolve_pdb_chain_name(pdb_chain_id):
	cmd = "grep {} {}".format(pdb_chain_id,pdb_name_resolution_table)
	ret = subprocess.check_output(cmd, shell=True).split("\n")
	if len(ret)!=1: return None # resolution failed
	name_resolved = ret.split()[0].replace("_","")
	return name_resolved


##########################################
def find_pdb_ids_of_similar_seqs(sequence, cutoff_pct, scratch, tmpname):
	os.chdir(scratch)
	queryfile = tmpname+".fasta"
	outf = open(queryfile,"w")
	outf.write(">{}\n{}\n".format(tmpname, sequence))
	outf.close()

	# outfmt 6 is tabular format; the column headers are given using outfmt 7:
	# query acc., subject acc., % identity, alignment length, mismatches,
	#      gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	target_pct_idtty = OrderedDict()
	cmd = "{} -db {} -query {} -outfmt 6".format(blastp, pdb_blast_db, queryfile)
	for line in subprocess.check_output(cmd, shell=True).split("\n"):
		print "\t\t", line
		field = line.split()
		if len(field)==0 or field[0] != tmpname:  continue # this is not the result line
		pct_identity = field[2]
		if float(pct_identity)<cutoff_pct: break
		target_id = field[1].replace("_","")
		if target_id[-1]==target_id[-2]: target_id=target_id[:-1]
		if len(target_id)>5:
			target_id = resolve_pdb_chain_name(field[1])
			if not target_id: continue # name resolving failed for one reason or another

		if not target_id in target_pct_idtty.keys(): target_pct_idtty[target_id]=pct_identity
	return target_pct_idtty
