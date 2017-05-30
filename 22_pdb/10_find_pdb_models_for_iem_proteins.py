#!/usr/bin/python

from integrator_utils.mysql import *


##########################################
def main():

	db, cursor = connect()

	qry = 'select gene_name from omim_genemaps where inborn_error_of_metabolism=1'
	ret = search_db(cursor, qry)
	for line in ret:
		print line


#########################################
if __name__ == '__main__':
	main()
