#!/usr/bin/python3
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from integrator_utils.python.mysql import *


def bed_format(data, start, end, exon_intron_label = None):
    columns = [data["chrom"],  start, end,  data["name"], data["name2"], data["strand"]]
    if exon_intron_label: columns.append(exon_intron_label)
    return "\t".join([str(c) for c in columns])


def csv_format(data, start, end, exon_intron_label):
    columns = [data["name"], data["name2"], data["strand"], start, end]
    if exon_intron_label: columns.append(exon_intron_label)
    return "\t".join([str(c) for c in columns])


def output_gene(data, outformatter, outf):
    leftmost  = data["exonStarts"].split(',')[0]
    rightmost = data["exonEnds"].split(',')[-2]  # ucsc exon starts/ends end with ','
    print(outformatter(data, leftmost, rightmost), file=outf)


def output_exons(data, outformatter, outf):

    starts = data["exonStarts"].split(',')[:-1]
    ends   = data["exonEnds"].split(',')[:-1]  # ucsc exon starts/ends end with ','

    if len(starts) != len(ends):
        print("start / end lenght mismatch\n", data)
        exit(1)

    exon_ct  = 0
    prev_end = None
    output = []
    for start, end in zip(starts, ends):
        if prev_end:
           output.append(outformatter(data, prev_end, start, f"intron_{exon_ct}"))
        exon_ct += 1
        output.append(outformatter(data, start, end, f"exon_{exon_ct}"))
        prev_end = end

    for line in output:
        print(line, file=outf)


#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    outformat  = "bed"
    outcontent = "exons"
    db     = connect_to_mysql(conf_file="/home/ivana/.ucsc_mysql_conf")
    if not db: exit(1)
    cursor = db.cursor()

    table = "refGene"
    cols_to_extract = ["chrom", "name", "name2", "strand", "exonStarts", "exonEnds"]
    qry  = f"select {','.join(cols_to_extract)} from {table} "
    outdir = "/storage/databases/ucsc/gene_regions"

    # for database in ['hg17', 'hg18', 'hg19']:
    outformatter = bed_format if outformat == "bed" else csv_format
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
            chrom = data["chrom"]
            if chrom not in outhandle:
                outhandle[chrom] = open(f"{outdir}/{database}/{chrom}.{outformat}", "w")

            if outcontent == "exons":
                output_exons(data, outformatter, outhandle[chrom])
            else:
                output_gene(data, outformatter, outhandle[chrom])

        for outf in outhandle.values(): outf.close()

    cursor.close()
    db.close()
    
    return True


#########################################
if __name__ == '__main__':
    main()


