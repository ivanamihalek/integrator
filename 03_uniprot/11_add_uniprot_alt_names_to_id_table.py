#! /usr/bin/python

from integrator_utils.mysql import *
from sys import argv
from re import match, search, sub
from socket import gethostname
development =  gethostname()=='pegasus'

####################################################
def find_gene_by_uniprot_id(cursor, uniprot_id):
	qry = 'select id,  symbol, synonyms, uniprot_ids from genes where uniprot_ids like "%%%s%%"' % uniprot_id
	ret = search_db(cursor,qry)
	if ret: return ret
	# otherwise, we do not have the uniprot_id properly stored in genes table
	# we are not going to fix this here
	# find gene by the gene_name(symbol)
	qry = "select  gene_name from uniprot_basic_infos where uniprot_id='%s'" % uniprot_id
	ret = search_db(cursor,qry)
	if not ret: return None
	# uniprot_id is primary key in uniprot_basic_infos so it should be unique
	# (i.e. there should not be more than one name here) - I am pretending that's not my problem here
	gene_name = ret[0][0]
	qry = 'select id,  symbol, synonyms, uniprot_ids from genes where symbol="%s"' % gene_name
	ret = search_db(cursor,qry)
	if not ret: return None
	# here however, if there is more than 1 return, I bail out
	if len(ret)!=1: return None
	# otherwise, I hope that's it
	return ret

####################################################
def store(cursor, uniprot_id, gene_name, uniprot_synonyms):
	ret = find_gene_by_uniprot_id(cursor, uniprot_id)
	if not ret: return
	for row in ret:
		new_alias_symbols = uniprot_synonyms
		[id, symbol, existing_synonyms, uniprot_ids] = row
		if gene_name!= symbol: continue # it should not happen
		if not uniprot_ids: continue # not sure how this happens, but it does
		if existing_synonyms and existing_synonyms != "":
			new_alias_symbols += existing_synonyms.replace(' ','').replace('\'','\'\'').split("|")
		qry  = "update genes set synonyms='%s' where id=%d" % ("|".join(list(set(new_alias_symbols))), id)
		search_db(cursor,qry,verbose=True)
	return

####################################################
def parse(entry):
	(uniprot_id, gene_name) = ("", "")
	synonyms= []
	if len(entry)==0: return [uniprot_id, gene_name, synonyms]
	if not search(r'\nOS   Homo sapiens', entry): return [uniprot_id, gene_name, synonyms]

	# the Name field may be split over several lines starting with GN
	# it also has xrefs in {} (to get rid of)
	name_lines = ""
	for line in  entry.split("\n"):
		if (line[:2] == "AC"):
			# it looks like the first iD is current,
			# and the rest are the older ids
			m = match(r'AC   (\w{6});',line)
			if not m: continue
			uniprot_id = m.group(1)
			pass
		elif (line[:2]=="GN"):
			name_lines += line

	name_lines = sub(r'\{.*?\}','', name_lines)
	# yes, there are cases without a name; see Q1T7F1
	if name_lines == "": return [uniprot_id, gene_name, synonyms]

	m = search(r'Name=([\w\-\s]+?);', name_lines)
	if not m: return [uniprot_id, gene_name, synonyms]
	gene_name = m.group(1).replace(' ','')

	if "SPAR" in synonyms:
		print entry
		exit()

	m = search(r'Synonyms=([\w\,\s\-]+);', name_lines)
	if m:  synonyms = m.group(1).replace(' ','').replace('\'','\'\'').split(",")

	return [uniprot_id, gene_name, synonyms]


###################################################
def main():

	filename = "/databases/uniprot/uniprot_sprot.dat"

	if development:
		db = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
	else:
		db = connect_to_mysql(user="blimps", passwd=(os.environ['BLIMPS_DATABASE_PASSWORD']))
	if not db: exit(1)
	cursor = db.cursor()
	qry = 'set autocommit=1'  # not sure why this has to be done explicitly - it should be the default
	search_db(cursor, qry, False)
	if development:
		switch_to_db(cursor, 'blimps_development')
	else:
		switch_to_db(cursor, 'blimps_production')

	inf = open (filename, "r")

	entry = ""
	for line in inf:
		if line[:2]=="//":
			[uniprot_id, gene_name, synonyms] = parse(entry)
			if len(synonyms)>0:
				store(cursor, uniprot_id, gene_name, synonyms)
			entry = ""
		else:
			entry += line

	parse(entry)
	inf.close()

	cursor.close()

#########################################
if __name__ == '__main__':
	main()
