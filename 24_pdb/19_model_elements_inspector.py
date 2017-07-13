#!/usr/bin/python
#
#
from integrator_utils.mysql     import *


##########################################
def newborn_screening_genes(cursor):
	nbs_genes = []
	switch_to_db(cursor,'monogenic_development')
	qry = "select name_short, omim_ids from diseases"
	for line in search_db(cursor, qry):
		name_short, omim_ids = line
		for omim_id in omim_ids.split(';'):
			qry  = "select approved_symbol, ensembl_gene_id from blimps_development.omim_genemaps "
			qry += "where mim_number='%s'" % omim_id
			ret = search_db(cursor,qry, verbose=True)
			approved_symbol, ensembl_gene_id = ret[0]
			qry  = "select uniprot_id, ec_number, cofactors from blimps_development.uniprot_basic_infos where ensembl_gene_id='%s'" % ensembl_gene_id
			ret = search_db(cursor,qry, verbose=True)
			uniprot_id, ec_number,uniprot_cofactors = ret[0]
			nbs_genes.append([approved_symbol, name_short, ensembl_gene_id, uniprot_id, ec_number, uniprot_cofactors])
	return nbs_genes

##########################################
def all_iem_related_genes(cursor):
	iem_genes = []
	qry  = "select mim_number, approved_symbol, ensembl_gene_id, phenotypes "
	qry += "from blimps_development.omim_genemaps "
	qry += "where inborn_error_of_metabolism=1"
	ret = search_db(cursor,qry, verbose=False)
	if not ret:
		print "no iems ?!"
		exit()
	for row in ret:
		[mim_number, approved_symbol, ensembl_gene_id, phenotypes] = row
		qry  = "select uniprot_id, ec_number, cofactors from blimps_development.uniprot_basic_infos where ensembl_gene_id='%s'" % ensembl_gene_id
		ret = search_db(cursor,qry, verbose=False)
		if not ret:
			print "no uniprot info found for", ensembl_gene_id
			continue
			#exit()
		uniprot_id, ec_number, cofactors = ret[0]
		#print mim_number, approved_symbol, ensembl_gene_id, uniprot_id, ec_number, phenotypes
		iem_genes.append([approved_symbol, phenotypes, ensembl_gene_id, uniprot_id, ec_number, cofactors])
	return iem_genes


##########################################
def main():
	db, cursor = connect()
	# first lets focus on the proteins from the newborn screening
	# genes = newborn_screening_genes(cursor)
	disease_descriptors = all_iem_related_genes(cursor)
	model_elements_found = 0
	total = 0
	enzymes = 0
	for dd in disease_descriptors:
		[gene_symbol, disease, ensembl_gene_id, uniprot_id, ec_number, uniprot_cofactors] = dd
		qry = "select * from monogenic_development.model_elements where gene_symbol='%s'" % gene_symbol
		ret = search_db(cursor,qry)
		if ec_number: enzymes += 1
		if ret: model_elements_found += 1
		total += 1

	print "total:", total, "    enzymes:", enzymes
	print "model_elements_found:", model_elements_found, " (%.2f)" % (float(model_elements_found)/float(enzymes))

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
