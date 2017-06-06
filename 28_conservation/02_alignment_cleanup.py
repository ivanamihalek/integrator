#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess

alignment_repository = "/home/ivana/monogenic/public/alignments"
afa2msf              = "/home/ivana/pypeworks/integrator/integrator_utils/afa2msf.pl"
remove_gap_only      = "/home/ivana/pypeworks/integrator/integrator_utils/remove_gap_only.pl"
patcher              = "/home/ivana/pypeworks/integrator/integrator_utils/patcher.pl"
muscle               = '/usr/bin/muscle'
##########################################
def main():

	for dependency in [alignment_repository, afa2msf, remove_gap_only, patcher,muscle]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()
	db, cursor = connect()
	qry = "select * from monogenic_development.diseases"
	ret = search_db(cursor, qry)
	for line in ret:
		[id, name_short, name_long, omim_ids, description, prim, sec] = line
		print "="*60
		print name_short
		for omim_id in omim_ids.split(";"):
			qry = "select approved_symbol, ensembl_gene_id from omim_genemaps where mim_number='%s'" % omim_id
			ret2 = search_db(cursor, qry)
			for line2 in ret2:
				approved_symbol = line2[0]
				ensembl_gene_id = line2[1]
				qry = "select uniprot_id from uniprot_basic_infos where gene_name='%s'" % approved_symbol
				ret3 = search_db(cursor, qry)
				if not ret3:
					print "uniprot id not found"
					exit()
				uniprot_id = ret3[0][0]
				qry = "select sequence from monogenic_development.uniprot_seqs where uniprot_id='%s'" % uniprot_id
				ret3 = search_db(cursor, qry)
				if not ret3 or len(ret3)==0:
					print "sequence not found"
					exit()
				sequence = ret3[0][0]
				print "\t", approved_symbol,ensembl_gene_id, uniprot_id
				uniprot_fa = "%s.fa"%uniprot_id
				outfa = open(uniprot_fa,"w")
				outfa.write(">%s\n%s\n"%(uniprot_id,sequence))
				outfa.close
				aln_name = "%s.afa"%ensembl_gene_id
				aln_path = alignment_repository+"/"+aln_name
				if not os.path.exists(aln_path):
					print aln_path, "not found"
					continue
				if os.stat(aln_path).st_size == 0:
					print aln_path, "empty"
					continue
				# get rid of Z's, tfm to msf, get rid of gaps only, patch
				# muscle and maft choke subprocess call - comething in off with python - move that piece to perl
				uniprot_patched_afa = "{}.patched.afa".format(uniprot_id)
				for cmd in ["sed 's/[\Z\*]/-/g' {} > tmp.afa".format(aln_path),
							"{} tmp.afa > tmp.msf".format(afa2msf),
							"{} tmp.msf > tmp2.msf".format(remove_gap_only),
							"{} tmp2.msf {} 0.62 > /dev/null".format(patcher, uniprot_patched_afa),
				            "rm tmp.* tmp2.*"]:
					subprocess.call(cmd, shell=True)

				os.rename(uniprot_fa, alignment_repository+"/"+uniprot_fa)
				os.rename(uniprot_patched_afa, alignment_repository+"/"+uniprot_patched_afa)


########################################
if __name__ == '__main__':
	main()
