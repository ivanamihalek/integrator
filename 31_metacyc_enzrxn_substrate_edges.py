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
def proces_reaction(cursor,  enzrxn_id, reaction_id):
	# I want to have an edge directly between enzrxn and left and right reaction compounds
	# 2 problems here:
	#   1) protons, x= oxygens and similar I do not want as nodes in the network
	#   2) there are classes of proteins and similar referred to here, which are not listed explicitly
	#
	# for 1 use histogram from the previous script
	# for 2: there is a class table in metacyc; still, not clear which proteins fore example begin to a class
	qry  = "select rxn_left, rxn_right from metacyc_reactions where unique_id like '%s'" % reaction_id
	rows = search_db(cursor, qry)
	print rows
	return

##########################################
def main():
	db, cursor = connect()

	enzrxn_ids = get_enzrxn_ids(db, cursor)
	for enzrxn_id in enzrxn_ids:
		qry  = "select reaction, cofactors, required_protein_complex, "
		qry += "regulated_by, alternative_substrates, alternative_cofactors "
		qry += "from metacyc_enzrxns where unique_id='%s'" % enzrxn_id
		rows = search_db(cursor, qry)
		if not rows or len(rows)==0: continue
		for row in rows:
			[reaction, cofactors, required_protein_complex, regulated_by,
			 alternative_substrates, alternative_cofactors] = row
			print '===============', enzrxn_id
			print reaction
			proces_reaction(cursor, enzrxn_id, reaction)

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
