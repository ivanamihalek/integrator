#!/usr/bin/python

# the core table in metacyc seems to be enzrxns
# it uses field | to link to

# enzyme                   | protein     --->  components | protein
# reaction                 | pathway     --->  in_pathway | pwy

from integrator_utils.mysql import *
import os

##########################################
def get_pathway_ids(db, cursor):
	qry  = "select distinct(to_id) from  metacyc_edges where to_field='pathway'"
	rows = search_db(cursor, qry)
	if not rows or len(rows)==0:
		print "no pathways found (?!)"
		cursor.close()
		db.close()
		exit()

	return [row[0] for row in rows]

##########################################
def store_pwy_link (cursor, sub_pwy_id, super_pwy_id):
	fixed_fields = {'from_field': 'pathway', 'from_id': sub_pwy_id,
					'to_field': 'pathway','to_id': super_pwy_id,
					'via': 'sub_pathway'}
	update_fields = {}
	store_or_update (cursor, 'metacyc_edges', fixed_fields, update_fields, verbose=True)
	return

##########################################
def walk_both_directions(cursor, pwy_id, visited):

	if pwy_id in visited: return

	visited.append(pwy_id)

	# let the database keep track if we have already seen this pair, just keep it in the order sub to super
	qry = "select in_pathway, sub_pathways, super_pathways "
	qry += "from metacyc_pathways where unique_id= '%s' " % pwy_id
	rows = search_db(cursor, qry)
	if not rows or len(rows)==0:
		return
	for row in rows:
		[in_pathway, sub_pathways, super_pathways] = row
		if sub_pathways!= '':
			for sub_pwy in sub_pathways.replace(' ','').split(";"):
				if sub_pwy =="": continue
				store_pwy_link (cursor, sub_pwy, pwy_id)
				walk_both_directions(cursor, sub_pwy, visited)
		supers = set()
		if sub_pathways != '':
			supers |= set(super_pathways.replace(' ','').split(";"))
		if in_pathway != '':
			supers |= set(in_pathway.replace(' ','').split(";"))
		supers -= set([""])
		for super_pwy in supers:
			store_pwy_link (cursor, pwy_id, super_pwy)
			walk_both_directions(cursor, super_pwy, visited)

	return

##########################################
def main():
	db, cursor = connect()

	pathway_ids = get_pathway_ids(db, cursor)
	visited = []
	for pwy_id in pathway_ids:
		walk_both_directions(cursor, pwy_id, visited)

	print len(visited)
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()

