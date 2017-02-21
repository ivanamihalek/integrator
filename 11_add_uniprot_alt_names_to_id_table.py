#! /usr/bin/python

from integrator_utils.mysql import *
from sys import argv
from re import match, search

####################################################
def find_gene_by_uniprot_id(cursor, uniprot_id):
	qry = 'select id,  symbol, alias_symbol, uniprot_ids from genes where uniprot_ids like "%%%s%%"' % uniprot_id
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
	qry = 'select id,  symbol, alias_symbol, uniprot_ids from genes where symbol="%s"' % gene_name
	ret = search_db(cursor,qry)
	if not ret: return None
	# here however, if there is more than 1 return, I bail out
	if len(ret)!=1: return None
	# otherwise, I hope that's it
	return ret

####################################################
def store(cursor, uniprot_id, gene_name, synonyms):
	print  uniprot_id, gene_name, synonyms
	ret = find_gene_by_uniprot_id(cursor, uniprot_id)
	if not ret: return
	for row in ret:
		new_alias_symbols = []
		[id, symbol, alias_symbol, uniprot_ids] = row
		if not uniprot_ids: continue # not sure how this happens, but it does
		# paranoid android
		if not uniprot_id in uniprot_ids.split('|'): continue
		if alias_symbol and alias_symbol != "": new_alias_symbols = alias_symbol.replace(' ','').split('|')
		for syn in synonyms:
			# make sure our 'alias' not actually the official symbol
			# also make sure we don't have it already
			if syn in [symbol]+new_alias_symbols: continue
			new_alias_symbols.append(syn.replace(' ','')) #not sure where these spaces are coming from
		new_alias_symbol = '|'.join(new_alias_symbols).replace(' ','')
		if alias_symbol == new_alias_symbol: continue
		qry  = "update genes set alias_symbol='%s' " % new_alias_symbol
		qry += "where id=%d" % id
		search_db(cursor,qry,verbose=True)
	return

####################################################
def parse(entry):
	(uniprot_id, gene_name) = ("", "")
	synonyms= []
	if len(entry)==0: return [uniprot_id, gene_name, synonyms]
	if not search(r'\nOS   Homo sapiens', entry): return [uniprot_id, gene_name, synonyms]


	for line in  entry.split("\n"):
		if (line[:2] == "AC"):
			# it looks like the first iD is current,
			# and the rest are the older ids
			m = match(r'AC   (\w{6});',line)
			if not m: continue
			uniprot_id = m.group(1)
			pass
		elif (line[:10]=="GN   Name="):
			m = search(r'Name=(\w+);*', line)
			gene_name = m.group(1)
			m = search(r'Synonyms=([\w\,\s]+);*', line)
			if not m: continue
			synonyms = m.group(1).replace(r'\s','').split(",")

	return [uniprot_id, gene_name, synonyms]



###################################################
def main():
	if len(argv) < 2:
		print  "Usage: %s  <file name>" % argv[0]
		exit(1)

	#db = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
	db = connect_to_mysql(user="blimps_production", passwd=(os.environ['BLIMPS_DATABASE_PASSWORD']))
	if not db: exit(1)
	cursor = db.cursor()
	qry = 'set autocommit=1'  # not sure why this has to be done explicitly - it should be the default
	search_db(cursor, qry, False)
	switch_to_db(cursor, 'blimps_development')

	filename = argv[1]
	inf = open (filename, "r")

	entry = ""
	for line in inf:
		if line[:2]=="ID":
			[uniprot_id, gene_name, synonyms] = parse(entry)
			if len(synonyms)>0: store(cursor, uniprot_id, gene_name, synonyms)
			entry = ""
		else:
			entry += line

	parse(entry)
	inf.close()

	cursor.close()

#########################################
if __name__ == '__main__':
	main()
