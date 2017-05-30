#!/usr/bin/python

from integrator_utils.mysql import *

# the origin list of metabolic disorders with OMIM number was found here
# http://www.ssiem.org/centralstore/resources/SSIEMClassificationIEM2011.pdf
# the omim numbers were extracted by some crude copying to textfile and grepping for all-digit strings
# (the resulting set of number, at the time of this writing,  was soleved in
#  omim_inherited_metabolic_disorders_list.txt, as a single column

# the following numbers seem to be misassigned (sporadic vs familial type of disease and similar)
replace = {}
replace["176090"] = "176100"
replace["194470"] = None
# HYPERZINCEMIA AND HYPERCALPROTECTINEMIA, 194470, seems to be related to albunmin,
# but OMIM does not have it in the genemap (?)
replace["201460"] = "201475" # moved by omim
replace["210100"] = None
replace["236130"] = None
replace["236795"] = None
replace["236900"] = None
replace["237400"] = None
replace["238340"] = None
replace["239400"] = None # removed from OMIM
replace["245130"] = None
replace["250951"] = None
replace["258920"] = "258900" # moved by omim
replace["612934"] = "614921" # moved by omim
replace["600737"] = "605820" # moved by omim
replace["307600"] = "602014" # moved by omim
replace["540000"] = None # mitochondrialy encoded
replace["545000"] = None # mitochondrialy encoded
replace["535000"] = None # mitochondrialy encoded
replace["530000"] = None
replace["520000"] = None # mitochondrialy encoded
replace["500009"] = None # mitochondrialy encoded
replace["500002"] = None # mitochondrialy encoded
replace["500001"] = None # mitochondrialy encoded
replace["551500"] = None # mitochondrialy encoded
replace["551000"] = None # mitochondrialy encoded
replace["557000"] = None # mitochondrialy encoded
replace["609055"] = "256730" #mapping by Ivana - the phenotype may be caused by other genes too
fix = {}

#  mim_number | gene_symbols  | gene_name  | approved_symbol |
# entrez_gene_id | ensembl_gene_id | comments | phenotypes | mouse_gene_symbol
fix["300371"] = '300371, \"ABCD1\", \"ATP-binding cassette (ABC) D1\", \"ABCD1\", '  # this fixes disease id "300100"
fix["300371"] += '215, \"ENSG00000101986\", \"\", '
fix["300371"] += '\"Adrenomyeloneuropathy, adult; Adrenoleukodystrophy (300100)\", \"Abcd1(MGI:1349215)\"'

fix["309060"] = '309060,  \"LAMP2 \", \"lysosome-associated membrane protein 2\", \"LAMP2 \", '
fix["309060"] += '\"3920\", \"ENSG00000005893\", \"\", '
fix["309060"] += '\"Danon disease (300257)\", \"Lamp2(MGI:96748)\"'

##########################################
def fix_omim_table(cursor):
	for mim_number in fix.keys():
		qry = 'select * from omim_genemaps where mim_number = %d ' % int(mim_number)
		ret = search_db(cursor, qry)
		if not ret:
			qry = "INSERT INTO omim_genemaps "
			qry +="(mim_number,gene_symbols,gene_name,approved_symbol,"
			qry += "entrez_gene_id,ensembl_gene_id,comments,phenotypes,mouse_gene_symbol)"
			qry +="VALUES (%s) " % fix[mim_number]
			ret = search_db(cursor, qry, verbose=True)



##########################################
def main():
	infile = open("/databases/omim/omim_inherited_metabolic_disorders_list.txt")

	db, cursor = connect()

	fix_omim_table(cursor)

	not_found = []
	total_found = 0
	total = 0
	for line in infile:
		total += 1
		omim_disease_id = line.rstrip()
		if replace.has_key(omim_disease_id):
			if not replace[omim_disease_id]: continue
			omim_disease_id = replace[omim_disease_id]
		qry = 'select * from omim_genemaps where phenotypes like "%%%s%%"' % omim_disease_id
		ret = search_db(cursor, qry)
		if not ret:
			# it looks like there is some mixup wiht disease and gene number
			qry = 'select * from omim_genemaps where mim_number = %d ' % int(omim_disease_id)
			ret = search_db(cursor, qry)
			if not ret:
				not_found.append(omim_disease_id)
				continue
		for line in ret:
			[id, mim_number , gene_symbols , gene_name , approved_symbol,
			entrez_gene_id , ensembl_gene_id , phenotypes , mouse_gene_symbol , inborn_error_of_metabolism] = line
			print id, mim_number, phenotypes
			qry = "update omim_genemaps set inborn_error_of_metabolism=1 where id=%d" %id
			search_db(cursor, qry)

		total_found += 1
	print "total diseases scanned:", total, "mapped to gene", total_found

	#for omim_disease_id in not_found:
	#	print omim_disease_id

#########################################
if __name__ == '__main__':
	main()
