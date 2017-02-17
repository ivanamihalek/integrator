#!/usr/bin/python

from   integrator_utils.mysql import *
import os

#########################################
def read_translation_table(filename):
    translation = {}
    inf = open(filename, "r")
    for line in inf:
        line = line.rstrip('\n')
        fields = line.split('\t')
        id = fields[0].replace(' ','')
        if not id or id=='': continue
        translation[id] = fields[1:]
    return translation

#########################################
def inversion_check(metacyc_gene2prot, metacyc_prot2gene, gene_id_translation):
    for gene_id, prots in metacyc_gene2prot.iteritems():
        proteins = prots.split(";")
        for protein in proteins:
            if gene_id==metacyc_prot2gene[protein]: continue
            print "name conflict:", gene_id_translation[gene_id]
            print gene_id, " ** ", protein, " ** ", metacyc_prot2gene[protein]


#########################################
def store_if(hash, key, value):
    if hash.has_key(key):
        # assuming keys are not super/substrings
        if value in hash[key]: return
        hash[key] += ";"+value
    else:
        hash[key] = value
    return

#########################################
def main():

    protein_id_translation = read_translation_table("/databases/meta_cyc/data/datfiles/protein-links.dat")
    gene_id_translation = read_translation_table("/databases/meta_cyc/data/datfiles/gene-links.dat")

    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    db     = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
    if not db: exit(1)
    cursor = db.cursor()
    qry = 'set autocommit=1' # not sure why this has to be done explicitly - it should be the default
    search_db(cursor,qry,False)

    # find proteins that are human
    switch_to_db(cursor, 'blimps_development')
    qry = 'select unique_id, gene  from metacyc_proteins where species="TAX-9606"'
    metacyc_prot2gene = {}
    for row in search_db(cursor,qry,False):
        [protein_id, gene] =  row
        if metacyc_prot2gene.has_key(protein_id):
            print "two genes for a protein? protein:", protein_id, "genes:", metacyc_prot2gene[protein_id], gene
            exit(1)
        else:
            metacyc_prot2gene[protein_id] = gene.replace(' ','')
    # get gene2prot, the same way
    qry = 'select g.unique_id, g.product from metacyc_genes as g, metacyc_proteins as p '
    qry += ' where g.product= p.unique_id and p.species="TAX-9606"'
    metacyc_gene2prot = {}
    for row in search_db(cursor,qry,False):
        [gene_id,  protein] =  row
        store_if(metacyc_gene2prot, gene_id, protein.replace(' ',''))

    # inversion check:
    inversion_check(metacyc_gene2prot, metacyc_prot2gene, gene_id_translation)

    missed_keys = {}
    metacyc_uniprot2protein = {}
    metacyc_uniprot2protein["O94925"] = "MONOMER-11425"
    metacyc_uniprot2protein["P00747"] = "MONOMER-14919;MONOMER-14920"
    metacyc_uniprot2protein["F8WCM5"] = "ENSG00000129965-MONOMER"
    metacyc_uniprot2protein["P01308"] = "MONOMER-16191"

    for protein_id in  metacyc_prot2gene.keys():
        if not protein_id_translation.has_key(protein_id):
            missed_keys[protein_id] = "not in protein-links"
            continue
        if protein_id_translation[protein_id][1]=='':
            gene_id = protein_id_translation[protein_id][0]
            if gene_id=='':
                gene_id = metacyc_prot2gene[protein_id]
            if gene_id=='':
                missed_keys[protein_id] = "no gene found"
                continue
            if not gene_id_translation.has_key(gene_id):
                missed_keys[protein_id] = "no uniprot found for gene"
                continue
            store_if(metacyc_uniprot2protein, gene_id_translation[gene_id][0], protein_id)
        else:
            store_if(metacyc_uniprot2protein, protein_id_translation[protein_id][1], protein_id)

    #for id, reason in missed_keys.iteritems():
    #    print id, reason

    print " metacyc_prot2gene:", len(metacyc_prot2gene)
    print " metacyc_gene2prot:", len(metacyc_gene2prot)
    print " uniprot ids mapped to metacyc protein:", len(metacyc_uniprot2protein)
    print
    # scan through gene table
    # use uniprot ids to find metacyc_protein and metacyc_genes_ids
    # store them in gene table (add the columns to the table if needed)
    qry = "select id, symbol, name, uniprot_ids from genes"
    for row in search_db(cursor, qry):
        [id, symbol, name, uniprot_ids] = row
        if not uniprot_ids: continue
        wrote_gene = False
        for upid in uniprot_ids.split(";"):
            if not  metacyc_uniprot2protein.has_key(upid): continue
            if not wrote_gene: print id, symbol, name, uniprot_ids
            wrote_gene = True
            if not metacyc_uniprot2protein[upid] or  metacyc_uniprot2protein[upid]=='': continue
            for mc_protein_id in metacyc_uniprot2protein[upid].split(";"):
                if  not metacyc_prot2gene[mc_protein_id] or metacyc_prot2gene[mc_protein_id]=='': continue
                for mc_gene_id in metacyc_prot2gene[mc_protein_id].split(";"):
                    print "\t", mc_protein_id, " ***  ", mc_gene_id
                    print "\t", protein_id_translation[mc_protein_id][:3]
                    if mc_gene_id !='': print "\t", gene_id_translation[mc_gene_id]
        #if wrote_gene: exit()

    return True


#########################################
if __name__ == '__main__':
    main()


###################################################
# the following were found by hand
#  MONOMER-15797  GAS glucuronate 2-sulfatase small subunit
# [chondroitin-sulfate]-2-O-sulfo-b-D-glucuronate + H2O
#                   -> [chondroitin]-b-D-glucuronate + sulfate + H+
#  I cannot find this enzyme in human - there was a single paper in 1985,
# and all links take me back to it
#  MONOMER-15798  GAS  glucuronate 2-sulfatase large subunit
#  MONOMER-14910 not found in metacyc online
#  MONOMER-19878 not found in metacyc online
#  MONOMER-16542 not found in metacyc online
#  MONOMER-19876 not found in metacyc online
#  MONOMER-19877 not found in metacyc online
#  G66-33951-MONOMER  might be zebra fish
#  MONOMER-13386 not found in metacyc online
#  MONOMER-16018 Streptomyces venezuelae
#  MONOMER-12921 not found in metacyc online
#  MONOMER-19881 not found in metacyc online
#  MONOMER-19880 not found in metacyc online
#  MONOMER-14907 not found in metacyc online

#  MONOMER-11425 GLS 	O94925    glutaminase kidney isoform large subunit,  L-glutamine + H2O = L-glutamate + NH3
#  MONOMER-14919 PLG    P00747   plasminogen, plasmin large subunit [a fibrin + H2O -> n a peptide
#  MONOMER-14920 PLG             plasmin small subunit
#                 Plasmin is released as a zymogen called plasminogen (PLG)
#                 from the liver into the factor IX systemic circulation and
#                 placed into the MD5+ that leads into the lungs.
#                 Plasminogen is the precursor to the fibronolytic protease
#                 plasmin as well as angiostatin, an angiogenesis/tumor metastasis inhibitor.
# plasminogen processing:
# http://www.sigmaaldrich.com/life-science/metabolomics/enzyme-explorer/analytical-enzymes/plasmin.html

#  ENSG00000129965-MONOMER  INS-IGF2  F8WCM5 insulin isoform 2
#  MONOMER-16191 INS  P01308  insulin B chain


