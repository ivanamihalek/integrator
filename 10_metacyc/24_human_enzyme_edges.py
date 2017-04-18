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
import os


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
	qry  = 'select unique_id, enzyme from  metacyc_enzrxns'
	rows = search_db(cursor, qry)
	for row in rows:
		[enzrxn_id, enzyme_unique_id] = row
		monomers = find_human_monomers(cursor, enzyme_unique_id)
		if len(monomers)==0: continue
		proteins[enzrxn_id] = list(set(monomers))

	print len(proteins.keys())
	for enzrxn, proteins in proteins.iteritems():
		for prot in proteins:
			fixed_fields= {'from_field' : 'enzrxn',
						   'from_id' : enzrxn,
						   'to_field': 'protein',
						   'to_id' :prot,
						   'via' : 'enzyme'}
			update_fields = {}
			store_or_update (cursor, 'metacyc_edges', fixed_fields, update_fields, verbose=True)
			#exit()

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()
