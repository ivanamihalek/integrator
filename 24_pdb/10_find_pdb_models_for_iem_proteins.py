#!/usr/bin/python

from integrator_utils.mysql import *
import os
import urllib

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
def parse_url(ur):
	# url of the form
	#https://swissmodel.expasy.org/repository/uniprot/P22304.pdb?from=34&to=550&template=5fql.1.A&provider=swissmodel
	[addr, fnm] = ur.split('?')
	uniprot_id = addr.split("/")[-1].replace('.pdb','')
	pdbname =  uniprot_id
	for field in fnm.split("&"):
		[k,v] = field.split("=")
		pdbname += "_" + v
	pdbname += ".pdb"
	return pdbname
##########################################
def main():
	swissmodel_dir = "/databases/swissmodel"
	swissmodel_meta_file = swissmodel_dir+"/index_human.csv"

	db, cursor = connect()

	pdbmodels = parse_meta(swissmodel_meta_file)

	#qry = 'select approved_symbol, ensembl_gene_id, phenotypes from omim_genemaps where inborn_error_of_metabolism=1'
	qry = "select approved_symbol, ensembl_gene_id, phenotypes from omim_genemaps where approved_symbol='MUT'"
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
		qry  = "select uniprot_id, canonical_aa_length, old_ids from uniprot_basic_infos where ensembl_gene_id like '%%%s%%'" % ensembl_gene_id
		ret2 = search_db(cursor, qry)
		proteins += len(ret2)
		if not ret2:
			uniprot_not_found.append(ensembl_gene_id)
			continue
		for line2 in ret2:
			[uniprot_id, canonical_aa_length, old_ids] = line2
			if not pdbmodels.has_key(uniprot_id):
				uniprot_id_found = None
				for old_id in old_ids.split(";"):  # looks like this is never really the problem
					if pdbmodels.has_key(old_id):
						uniprot_id_found = old_id
						break
				if not uniprot_id_found:
					no_models.append(uniprot_id)
					continue
				uniprot_id = uniprot_id_found
			# make folder if it does not exist
			path = "/".join([swissmodel_dir,approved_symbol[0], approved_symbol])
			if not os.path.exists(path): os.makedirs(path)
			for model in pdbmodels[uniprot_id]:
				[coordinate_id, provider, start, end, template, qmean, qmean_norm,url] = model
				fraction = float(int(end)-int(start)+1)/int(canonical_aa_length)>0.75
				pdbname = parse_url(url)
				filepath = path+"/"+pdbname
				if not os.path.exists(filepath) or os.stat(filepath).st_size == 0: urllib.urlretrieve (url, filepath)

	print "proteins:", proteins
	if len(uniprot_not_found)>0:
		print 'uniprot entries not found:'
		for ensembl_gene_id in uniprot_not_found:
			print ensembl_gene_id

	if len(no_models)>0:
		print len(no_models),'pdb models not found:'
		for uniprot_id in no_models:
			print uniprot_id
	print "structure model found:", found

#########################################
if __name__ == '__main__':
	main()
