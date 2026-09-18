#! /usr/bin/python
# yet another patch ...

from integrator_utils.mysql import *
import os, subprocess
#from socket import gethostname
#development =  ( gethostname() in ['bison','pegasus'] )

###################################################
def main():

    assembly = 'hg19'
    ucsc_gene_regions_dir = "/databases/ucsc/canonical_gene_coords/%s" % assembly
    if not os.path.exists(ucsc_gene_regions_dir):
        print ucsc_gene_regions_dir, "not found"
        exit()

    db, cursor = connect()
    switch_to_db(cursor,'monogenic_development')
    qry = "select uniprot_id, ensembl_transcript_id from uniprot_seqs where hg19_exon_starts is null"
    ret = search_db(cursor,qry)
    if not ret:
        print "all uniprot_seqs appear to have hg19_exon_starts column filled"
        exit()
    if 'error' in ret[0][0].lower():
        search_db(cursor,qry, verbose=True)
        exit()
    for line in ret:
        [uniprot_id, ensembl_transcript_id] = line
        qry  = "select gene_name from blimps_development.uniprot_basic_infos "
        qry += "where uniprot_id='%s' " % uniprot_id
        ret = search_db(cursor,qry)
        if not ret or 'error' in ret[0][0].lower():
            print "gene name not found for '%s' " % uniprot_id
            exit()
        gene_name = ret[0][0]
        print uniprot_id, ensembl_transcript_id, gene_name
        cmd = "grep  %s %s/*" % (gene_name, ucsc_gene_regions_dir)
        ret =  subprocess.check_output(cmd, shell=True).rstrip()
        if not ret or len(ret)==0:
            print "no entry for %s found in %s " % (gene_name, ucsc_gene_regions_dir)
            continue
        lines = []
        for line in ret.split("\n"):
            fields = line.split("\t")
            [infile, gene_names] = fields[0].split(":")
            if not gene_name in gene_names.replace(" ","").split(","): continue
            lines.append(fields[1:])
        if len(lines)==0:
            print "no entry for %s found in %s " % (gene_name, ucsc_gene_regions_dir)
            continue
        if len(lines)==2:
            print "more than one entry found for %s found in %s " % (gene_name, ucsc_gene_regions_dir)
            continue
        # we assume certain format in the file name, containting the chromosome number: e.g. chr18.csv
        chromosome = infile.split("/")[-1].replace("chr","").replace(".csv","")
        [ucsc_id1, ucsc_id2, strand, txStart, txEnd, exonStarts, exonEnds] = lines[0]
        print chromosome, strand, txStart, txEnd, exonStarts, exonEnds
        if not exonStarts or len(exonStarts)==0:
            print "no entry for exonStarts (?)"
            exit(1)
        qry = "update uniprot_seqs "
        qry += "set "
        qry += "chrom='%s', strand='%s', " % ( chromosome, strand)
        qry += "hg19_exon_starts='%s', hg19_exon_ends='%s', " % (exonStarts, exonEnds)
        qry += "hg19_cds_start=%s, hg19_cds_end=%s " % (txStart, txEnd)
        qry += "where uniprot_id='%s' " % uniprot_id
        search_db(cursor,qry, verbose=True)

    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
