#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from   integrator_utils.mysql import *

#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    db     = connect_to_mysql(conf_file="/home/ivana/.ucsc_mysql_conf")
    if not db: exit(1)
    cursor = db.cursor()

    for database in ['hg17', 'hg18', 'hg19']:
        print "downloading from", database
        switch_to_db(cursor, database) # mouse build name
        # no Y chromosome, we are looking at uterus tissue
        outf = open("/databases/ucsc/%s_refgene.csv" % database, "w")
        qry  = "select  * from refGene  " 
        rows = search_db(cursor,qry)        
        for row in rows:
            print  >>outf,  "\t".join( [ str(r) for r in row] )
        outf.close()
        
    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


