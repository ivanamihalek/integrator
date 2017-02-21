#! /usr/bin/python

from integrator_utils.mysql import *

###################################################
def main():

	db = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
	if not db: exit(1)
	cursor = db.cursor()
	qry = 'set autocommit=1'  # not sure why this has to be done explicitly - it should be the default
	search_db(cursor, qry, False)
	switch_to_db(cursor, 'blimps_development')

	qry = "select id, alias_symbol from genes where not alias_symbol is  null"
	for row in search_db(cursor,qry):
		[id, alias_symbol] = row
		if not alias_symbol or alias_symbol == "": continue
		# for some reason replace(r'\s','') is not replacing space. I hate python.
		new_alias_symbol = '|'.join(list(set(alias_symbol.replace(' ','').upper().split('|'))))
		if alias_symbol == new_alias_symbol: continue
		qry  = "update genes set alias_symbol='%s' " % new_alias_symbol
		qry += "where id=%d" % id
		search_db(cursor,qry,verbose=True)
	cursor.close()

#########################################
if __name__ == '__main__':
	main()
