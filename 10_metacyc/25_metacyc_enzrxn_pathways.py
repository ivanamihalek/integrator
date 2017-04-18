#!/usr/bin/python

# the core table in metacyc seems to be enzrxns
# it uses field | to link to

# enzyme                   | protein     --->  components | protein
# reaction                 | pathway     --->  in_pathway | pwy

from integrator_utils.mysql import *
import os

##########################################
def get_enzrxn_ids(db, cursor):
	qry  = "select distinct(from_id) from  metacyc_edges where from_field='enzrxn'"
	rows = search_db(cursor, qry)
	if not rows or len(rows)==0:
		print "no human enrxns found (?!)"
		cursor.close()
		db.close()
		exit()

	return [row[0] for row in rows]

##########################################
def process_pathways(cursor, enzrxn_id, pathways):

	for pwy in pathways.replace(' ','').split(";"):
		fixed_fields= {'from_field' : 'enzrxn',
						'from_id' : enzrxn_id,
						'to_field': 'pathway',
						'to_id' : pwy,
						'via' : 'reaction'}
		update_fields = {}
		store_or_update (cursor, 'metacyc_edges', fixed_fields, update_fields, verbose=True)
	return

##########################################
def main():
	db, cursor = connect()

	reaction_substrates = []
	enzrxn_ids = get_enzrxn_ids(db, cursor)
	no_pathway = 0
	for enzrxn_id in enzrxn_ids:
		# also process alternative substrates here
		qry  =  "select r.in_pathway from metacyc_enzrxns as e, metacyc_reactions as r "
		qry  += "where e.unique_id='%s' and e.reaction = r.unique_id" % enzrxn_id
		rows = search_db(cursor, qry)
		if not rows or len(rows)==0 or len(rows[0][0].replace(' ',''))==0:
			no_pathway += 1
			continue
		if len(rows)>1:
			print 'caveat: did not expect more than a single reacton for enzrxn', enzrxn_id
			exit(1)

		process_pathways(cursor, enzrxn_id, rows[0][0])

	print "tot enzrxns:",  len(enzrxn_ids), " no pathay:",  no_pathway
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()

