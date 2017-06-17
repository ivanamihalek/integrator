#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess

pdb_dir        = "/databases/pdb"
swissmodel_dir = "/databases/swissmodel"
scratch_dir    = "/home/ivana/scratch"
scwrl          = "/usr/local/bin/scwrl4/Scwrl4"
pdb2seq       = "/home/ivana/pypeworks/integrator/integrator_utils/pdb2seq.pl"
swich_bb_identity = "/home/ivana/pypeworks/integrator/integrator_utils/switch_backbone_identity.pl"

cwd    = os.getcwd()
#########################################
def main():

	gene_symbol = "BTD"
	template_structure = "4cyg.pdb"
	template_full_path = "{}/{}".format(pdb_dir, template_structure)

	for dependency in [pdb_dir, swissmodel_dir, scratch_dir, scwrl, pdb2seq, swich_bb_identity, template_full_path]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	db, cursor = connect()
	switch_to_db(cursor, "monogenic_development")
	qry  = "select uniprot_seqs.sequence from uniprot_seqs, blimps_development.uniprot_basic_infos "
	qry += "where uniprot_seqs.uniprot_id=blimps_development.uniprot_basic_infos.uniprot_id "
	qry += "and blimps_development.uniprot_basic_infos.gene_name='%s'" % gene_symbol
	ret = search_db(cursor, qry)
	if not ret:
		print "seqeunce not found for", gene_symbol
		print "qry was:\n", qry

	target_seq = ret[0][0]
	cmd = "{} {}".format(pdb2seq, template_full_path)
	ret = subprocess.check_output(cmd, shell=True)
	if not ret:
		print "sequence not found for", template_full_path
		print "cmd was:\n", cmd
	template_seq = ret.rstrip()

	scratch = scratch_dir + "/" + gene_symbol + "_model"
	shutil.rmtree(scratch,ignore_errors=True)
	os.makedirs(scratch)
	os.chdir(scratch)
	outf = open("almt.fa","w")
	outf.write(">target\n%s\n"%target_seq)
	outf.write(">template\n%s\n"%template_seq)
	outf.close()
	#shutil.rmtree(scratch,ignore_errors=True)

	cursor.close()
	db.close()

 ########################################
if __name__ == '__main__':
	main()
