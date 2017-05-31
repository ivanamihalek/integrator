#! /usr/bin/python

from integrator_utils.mysql import *
from socket import gethostname
development =  gethostname()=='pegasus'

###################################################
def main():

	if development:
		db = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
	else:
		db = connect_to_mysql(user="blimps", passwd=(os.environ['BLIMPS_DATABASE_PASSWORD']))
	if not db: exit(1)
	cursor = db.cursor()
	qry = 'set autocommit=1'  # not sure why this has to be done explicitly - it should be the default
	search_db(cursor, qry, False)
	if development:
		switch_to_db(cursor, 'blimps_development')
	else:
		switch_to_db(cursor, 'blimps_production')

	qry = "select id, synonyms from genes where synonyms is not null"
	for row in search_db(cursor,qry):
		[id, alias_symbol] = row
		if not alias_symbol or alias_symbol == "": continue
		# for some reason replace(r'\s','') is not replacing space. I hate python.
		# to replace regexp, need to use re.sub - how hard is it to warn the user about such crap
		new_alias_symbol = '|'.join(list(set(alias_symbol.replace(' ','').replace("\'","\'\'").upper().split('|'))))
		if alias_symbol == new_alias_symbol: continue
		qry  = "update genes set synonyms='%s' " % new_alias_symbol
		qry += "where id=%d" % id
		search_db(cursor,qry,verbose=True)
	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
