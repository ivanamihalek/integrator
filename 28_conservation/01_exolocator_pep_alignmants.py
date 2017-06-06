#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess

alignment_repository = "/home/ivana/monogenic/public/alignments"

##########################################
def main():

	for dependency in [alignment_repository]:
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
				ensembl_gene_id  = line2[1]
				print "\t", approved_symbol,ensembl_gene_id
				aln_name =  "%s.afa"%ensembl_gene_id
				cmd = "wget  http://exolocator.bii.a-star.edu.sg/Best_MSA/pep/{}".format(aln_name)
				subprocess.call(cmd, shell=True)
				os.rename(aln_name, alignment_repository+"/"+aln_name)

########################################
if __name__ == '__main__':
	main()
