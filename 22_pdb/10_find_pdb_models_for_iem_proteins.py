#!/usr/bin/python

from integrator_utils.mysql import *

def parse_meta(swissmetafile):
	infile = open(swissmetafile,"r")
	models = {}
	for line in infile:
		if line[0]=='#': continue
		line = line.rstrip()
		[uniprot_id, coordinate_id, provider, start, end, template, qmean, qmean_norm,url] = line.split('\t')
		if not models.has_key(uniprot_id): models[uniprot_id] = []
		models[uniprot_id].append([coordinate_id, provider, start, end, template, qmean, qmean_norm,url])
	return models

##########################################
def main():

	swissmodel_meta_file = "/databases/swissmodel/index_human.csv"

	db, cursor = connect()

	pdbmodels = parse_meta(swissmodel_meta_file)

	qry = 'select approved_symbol, ensembl_gene_id, phenotypes from omim_genemaps where inborn_error_of_metabolism=1'
	ret = search_db(cursor, qry)
	print 'genes: ', len(ret)
	uniprot_not_found = []
	no_models = []
	found = 0
	proteins = 0
	for line in ret:
		[approved_symbol, ensembl_gene_id, phenotypes] = line
		if not ensembl_gene_id or ensembl_gene_id=="":
			print "no ensembl_gene_id: ", approved_symbol
			print phenotypes
			exit()
		print approved_symbol, ensembl_gene_id, phenotypes
		qry  = "select uniprot_id from uniprot_basic_infos where ensembl_gene_id like '%%%s%%'" % ensembl_gene_id
		ret2 = search_db(cursor, qry)
		proteins += len(ret2)
		if not ret2:
			uniprot_not_found.append(ensembl_gene_id)
			continue
		for line2 in ret2:
			[uniprot_id] = line2
			if not pdbmodels.has_key(uniprot_id):
				no_models.append(uniprot_id)
				print "no structure"
				continue

			print "\t", uniprot_id
			found += 1
			print pdbmodels[uniprot_id]

		exit()
	print "proteins:", proteins
	if len(uniprot_not_found)>0:
		print 'uniprot entries not found:'
		for ensembl_gene_id in uniprot_not_found:
			print ensembl_gene_id

	if False and len(no_models)>0:
		print 'pdb models not found:'
		for uniprot_id in no_models:
			print uniprot_id
	print "structure model found:", found

#########################################
if __name__ == '__main__':
	main()
