#!/usr/bin/python

# the core table in metacyc seems to be enzrxns
# it uses field | to link to

# enzyme                   | protein     --->  components | protein
# reaction                 | reaction    --->  rxn_left and rxn_right | compound
# cofactors                | compound
# required_protein_complex | protein (?)
# regulated_by             | regulation
# alternative_substrates   | compound
# alternative_cofactors    | compound


from integrator_utils.mysql import *
import os


##########################################
def get_enzrxn_ids(db, cursor):
	qry  = "select from_id from  metacyc_edges where from_field='enzrxn'"
	rows = search_db(cursor, qry)
	if not rows or len(rows)==0:
		print "no human enrxns found (?!)"
		cursor.close()
		db.close()
		exit()

	return [row[0] for row in rows]

##########################################
def proces_cofactors (cursor, enzrxn_id, cofactors):
	# again filter for too common
	print enzrxn_id, cofactors
	return

##########################################
def main():
	db, cursor = connect()

	enzrxn_ids = get_enzrxn_ids(db, cursor)
	for enzrxn_id in enzrxn_ids:
		qry  = "select cofactors, required_protein_complex, "
		qry += "regulated_by, alternative_substrates, alternative_cofactors "
		qry += "from metacyc_enzrxns where unique_id='%s'" % enzrxn_id
		rows = search_db(cursor, qry)
		if not rows or len(rows)==0:  continue
		for row in rows:
			[cofactors, required_protein_complex, regulated_by,
			 alternative_substrates, alternative_cofactors] = row
			if not cofactors or cofactors=='': continue
			print '===============', enzrxn_id
			print cofactors
			proces_cofactors (cursor, enzrxn_id, cofactors)

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
