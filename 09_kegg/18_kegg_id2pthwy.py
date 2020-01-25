#!/usr/bin/python3

from  integrator_utils.python.mysql import *
import os


def kegg_gene_ids_loaded(cursor):
	return hard_landing_search(cursor, "select count(*) from identifier_maps.kegg_human")[0][0] >10


def load_kegg_gene_ids(cursor):
	uni2kegg = "/storage/databases/uniprot/uniprot2keg.tsv"
	if not os.path.exists(uni2kegg)  or os.path.getsize(uni2kegg)==0:
		print(uni2kegg, "not found or is 0; run 17_uniprot_parser.pl first.")
		exit()
	inf = open(uni2kegg, "r")
	for line in inf:
		infields =  line.strip().split()
		if len(infields)<2: continue
		if len(infields) == 2:
			[uni, kegg] = infields
			old_uni = None
		else:
			[uni, kegg, old_uni] = infields
		if not kegg or len(kegg)==0 or not "hsa:" in kegg: continue
		fixed_fields = {"id":kegg.strip().replace("hsa:","")}
		update_fields = {"uniprot":uni}
		if old_uni:
			# dunno what kind of junk this is if there are more than 20 ids
			old_uni = ";".join(old_uni.split(";")[:20])
			update_fields["uniprot_old"]=old_uni
		store_or_update(cursor, 'kegg_human', fixed_fields, update_fields)

	inf.close()
	return


def load_kegg_pathways(cursor):
	kegg2kegg = "/storage/databases/kegg/pathway2gene.tsv"
	if not os.path.exists(kegg2kegg)  or os.path.getsize(kegg2kegg)==0:
		print(kegg2kegg, "not found; you can download it from http://rest.kegg.jp/link/hsa/pathway")
		exit()

	inf = open(kegg2kegg, "r")
	gene2path = {}
	path2gene = {}
	for line in inf:
		[pthwy,kegg] =  line.strip().split()
		pthwy = pthwy.replace("path:hsa", "")
		kegg = kegg.replace("hsa:", "")
		if not kegg in gene2path:
			gene2path[kegg]  = []
		gene2path[kegg].append(pthwy)
		if not pthwy in path2gene:
			path2gene[pthwy] = []
		path2gene[pthwy].append(kegg)

	print("number of genes mapped to pathway: ", len(gene2path))

	for kegg, pthwys in gene2path.items():
		pthwy_string = ";".join(sorted(pthwys, key=lambda x:len(path2gene[x])))
		# print(kegg, pthwy_string)
		# for pthwy  in pthwys:
		# 	print ("\t", pthwy, len(path2gene[pthwy]))
		fixed_fields = {"id":kegg.strip()}
		update_fields = {"kegg_pathways":pthwy_string}
		store_or_update(cursor, 'kegg_human', fixed_fields, update_fields)

	inf.close()
	return


def load_pathway_names(cursor):
	pthwy_names = "/storage/databases/kegg/pathway_names.txt"
	if not os.path.exists(pthwy_names)  or os.path.getsize(pthwy_names)==0:
		print(pthwy_names, "not found; you can download it from https://www.genome.jp/kegg-bin/get_htext")
		exit()
	inf = open(pthwy_names, "r")
	for line in inf:
		fields = line.strip().split()
		if len(fields)<2: continue
		if fields[0][0] != '0': continue
		print(fields[0], " ".join(fields[1:]))
		fixed_fields = {"id": fields[0]}
		update_fields = {"name":" ".join(fields[1:])}
		store_or_update(cursor, 'kegg_pathway_name', fixed_fields, update_fields)

	inf.close()
	return



##########################################
def main():
	mysql_conf_file = "/home/ivana/.tcga_conf"

	for dep in [mysql_conf_file]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	db = connect_to_mysql(conf_file=mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	#  get KEGG id: http://rest.kegg.jp/conv/genes/uniprot:P02763
	#  get pathways related to gene id:  http://rest.kegg.jp/link/pathway/hsa:10993
	#  get all genes related to pathway id: http://rest.kegg.jp/link/hsa/pathway:hsa01200
	#  link to  gene info https://www.kegg.jp/dbget-bin/www_bget?hsa:10993

	#
	switch_to_db(cursor, 'identifier_maps')
	if not kegg_gene_ids_loaded(cursor):
		print("kegg human ids do not seem to have been loaded")
		load_kegg_gene_ids(cursor)

	load_kegg_pathways(cursor)

	load_pathway_names(cursor)

	db.close()
	cursor.close()


#########################################
if __name__ == '__main__':
	main()
