#!/usr/bin/python
from integrator_utils.mysql import *
import os

##########################################
def compile_model(path,protein):
    print "compiling model for", protein

    exit()

##########################################
def main():
    swissmodel_dir = "/databases/swissmodel"
    #for subdir in list(map(chr, range(65,91))):
    #    path = "/". join([swissmodel_dir,subdir])
    #    if not os.path.exists(path): continue
    #    for protein in next(os.walk(path))[1]:
    #         compile_model(path,protein)
    db, cursor = connect()
    #switch_to_db(cursor, "monogenic_development")
    qry = "select * from monogenic_development.diseases"
    ret = search_db(cursor, qry)
    for line in ret:
        [id, name_short, name_long, omim_ids, description] = line
        print name_short
        for omim_id in omim_ids.split(";"):
            qry = "select approved_symbol from omim_genemaps where mim_number='%s'" % omim_id
            ret2 = search_db(cursor, qry)
            for line2 in ret2:
                protein = line2[0]
                print "\t", protein
                path = "/". join([swissmodel_dir,protein[0],protein])
                if not os.path.exists(path):
                    print path, "not found"
                    continue
                for model in next(os.walk(path))[2]:
                    print model
    return

#########################################
if __name__ == '__main__':
    main()
