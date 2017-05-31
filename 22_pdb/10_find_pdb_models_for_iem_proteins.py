#!/usr/bin/python

from integrator_utils.mysql import *


##########################################
def main():

	db, cursor = connect()

	qry = 'select approved_symbol, ensembl_gene_id, phenotypes from omim_genemaps where inborn_error_of_metabolism=1'
	ret = search_db(cursor, qry)
	uniprot_not_found = []
	for line in ret:
		[approved_symbol, ensembl_gene_id, phenotypes] = line
		print approved_symbol, ensembl_gene_id, phenotypes
		qry = "select uniprot_id from uniprot_basic_infos where ensembl_gene_id='%s'" % ensembl_gene_id
		ret2 = search_db(cursor, qry)
		if not ret2:
			uniprot_not_found.append(ensembl_gene_id)
			continue
		for line2 in ret2:
			[uniprot_id] = line2
			print "\t", uniprot_id


	for ensembl_gene_id in uniprot_not_found:
		print ensembl_gene_id
#########################################
if __name__ == '__main__':
	main()
