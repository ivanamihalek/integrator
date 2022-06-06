#!/usr/bin/python3
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from integrator_utils.python.mysql import *

#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    outformat = "bed"
    db     = connect_to_mysql(conf_file="/home/ivana/.ucsc_mysql_conf")
    if not db: exit(1)
    cursor = db.cursor()

    table = "refGene"
    cols_to_extract = ["chrom", "name", "name2", "strand", "exonStarts", "exonEnds"]
    qry  = f"select {','.join(cols_to_extract)} from {table} "
    outdir = "/storage/databases/ucsc/gene_regions"

    # for database in ['hg17', 'hg18', 'hg19']:
    for database in ['hg19']:
        print("downloading from", database)
        if not os.path.exists(f"{outdir}/{database}"): os.mkdir(f"{outdir}/{database}")
        colnames = get_column_names(cursor, database, "refGene")
        switch_to_db(cursor, database)  # mouse build name
        for col in cols_to_extract:
            if col in colnames: continue
            print(f"{col} not found among columns in {database}.{table}")
            exit(1)
        outhandle = {}
        rows = search_db(cursor, qry)
        for row in rows:
            data = dict(zip(cols_to_extract, [c.decode("utf-8") if type(c) == bytes else c for c in row]))
            leftmost   = data["exonStarts"].split(',')[0]
            rightmost = data["exonEnds"].split(',')[-2]  # ucsc exon starts/ends end with ','
            chrom = data["chrom"]
            if outformat == "bed":
                # ucsc coords are 0-based half open, as is the bed standard
                output = [chrom, leftmost, rightmost, data["name"], data["name2"], data["strand"]]
            else:
                output = [data["name"], data["name2"], data["strand"], leftmost, rightmost]

            if chrom not in outhandle:
                outhandle[data["chrom"]] = open(f"{outdir}/{database}/{chrom}.{outformat}", "w")
            print("\t".join([str(r) for r in output]), file=outhandle[chrom])

        for outf in outhandle.values(): outf.close()

    cursor.close()
    db.close()
    
    return True


#########################################
if __name__ == '__main__':
    main()


