#! /usr/bin/python

from integrator_utils.mysql import *
from sys import argv
from re import match, search, sub
from socket import gethostname
development =  gethostname()=='pegasus'

###################################################
def main():

    # read synonyms
    filename = "/databases/hgnc/aliases_and_prev_names.csv"
    inf = open(filename, "r")
    synonyms = {}
    for line in inf:
        onestr = ','.join(line.rstrip().split("\t"))
        fields = [x.replace(' ','') for x in onestr.split(',') if x !=""]
        symbol = fields[0]
        syns = "|".join(list(set(fields[1:])))
        synonyms[symbol] = syns
    inf.close

    # close
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

    qry = "select id, symbol from genes"
    for row in  search_db(cursor, qry, False):
        [id, symbol] = row
        if not row or row=="": continue
        if not synonyms.has_key(symbol) or not synonyms[symbol] or synonyms[symbol]=="": continue
        #print id, symbol, synonyms[symbol].replace('\'', '\'\'')
        qry  = "update genes set synonyms='%s' where id=%d " %  (synonyms[symbol].replace('\'', '\'\''), id)
        search_db(cursor, qry, True)


#########################################
if __name__ == '__main__':
	main()
