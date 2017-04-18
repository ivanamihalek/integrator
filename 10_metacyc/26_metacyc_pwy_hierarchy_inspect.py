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
def process_pwy_hierarchy(cursor, path_id):

	return

##########################################
def main():
	db, cursor = connect()

	pathway_ids = get_pathway_ids(db, cursor)

	no_further = 0
	has_sub    = 0
	in_and_super_the_same = 0
	not_the_same = 0
	for pwy_id in pathway_ids:
		qry = "select in_pathway, sub_pathways, super_pathways "
		qry += "from metacyc_pathways where unique_id= '%s' " % pwy_id
		rows = search_db(cursor, qry)
		if not rows or len(rows)==0:
			no_further += 1
			continue
		for row in rows:
			[in_pathway, sub_pathways, super_pathways] = row
			if in_pathway=='' and sub_pathways=='' and  super_pathways=='':
				no_further += 1
				continue
			if  sub_pathways!='':
				has_sub += 1
			# it looks like in and super are either the same thing, or one is subset of the other,
			# without much rhyme or reason - so just collect them into set
			if in_pathway==super_pathways or in_pathway in super_pathways or  super_pathways in in_pathway:
				in_and_super_the_same += 1
			else:
				not_the_same +=1
				print in_pathway, " *** ", super_pathways
	print "tot pathways:",  len(pathway_ids)
	print " no further :",  no_further
	print " has sub :",  has_sub
	print " in_and_super_the_same: ", in_and_super_the_same
	print " not the same:",  not_the_same

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()

