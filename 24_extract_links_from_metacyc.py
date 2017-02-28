#!/usr/bin/python

# the core table in metacyc seems to be enzrxns
# it uses field | to link to

# enzyme                   | protein     --->  components | protein
# reaction                 | reaction    --->  rxn_left and rxn_right | compound
# cofactor                 | compound
# required_protein_complex | protein (?)
# regulated_by             | regulation
# alternative_substrates   | compound
# alternative_cofactors    | compound


from integrator_utils.mysql import *
from socket import gethostname
import os


##########################################
def connect():
	development = gethostname()=='pegasus'
	if development:
		db = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
	else:
		db = connect_to_mysql(user="blimps", passwd=(os.environ['BLIMPS_DATABASE_PASSWORD']))
	if not db: exit(1)
	cursor = db.cursor()
	qry = 'set autocommit=1' # not sure why this has to be done explicitly - it should be the default
	search_db(cursor,qry,False)
	if development:
		switch_to_db(cursor, 'blimps_development')
	else:
		switch_to_db(cursor, 'blimps_production')
	return db, cursor

##########################################
def find_human_monomers(cursor, protein_unique_id):
	monomers = []
	qry = "select components,unmodified_form,species from metacyc_proteins where unique_id='%s'" % protein_unique_id
	rows =  search_db(cursor, qry)
	if not rows: return monomers
	[components,unmodified_form,species] = rows[0] # there shouldn't be two rows if the id is unique, should there
	if len(components)>0:
		for uniq_id in components.replace(' ','').split(";"):
			monomers += find_human_monomers(cursor, uniq_id)
	elif len(unmodified_form)>0:
		monomers += find_human_monomers(cursor, unmodified_form)
	elif species=='TAX-9606':
		monomers = [protein_unique_id]

	return monomers

##########################################
def main():
	db, cursor = connect()

	proteins = {}
	qry  = 'select id, enzyme from  metacyc_enzrxns'
	rows = search_db(cursor, qry)
	for row in rows:
		[enzrxn_id, enzyme_unique_id] = row
		monomers = find_human_monomers(cursor, enzyme_unique_id)
		if len(monomers)==0: continue
		print  enzrxn_id, list(set(monomers))
		proteins[enzrxn_id] = list(set(monomers))

	print len(proteins.keys())
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()

