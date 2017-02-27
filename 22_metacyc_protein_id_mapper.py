#!/usr/bin/python

# mysql> alter table metacyc_genes add column gene_id int(11);
# HS14864, protein  HS14864-MONOMER maps to deprecated identifier in Ensembl:
# this should be
# alpha-1,6-mannosylglycoprotein 6-beta-N-acetylglucosaminyltransferase V
# uniprot says that it is Q09328, MGAT5,
# forwhich we already hava a much healthier entry HS07793


from integrator_utils.mysql import *
from socket import gethostname
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
def mc_gene2prot_mapping(cursor):
    # get gene2prot, the same way
    qry = 'select g.unique_id, g.product from metacyc_genes as g, metacyc_proteins as p '
    qry += ' where g.product= p.unique_id and p.species="TAX-9606"'
    metacyc_gene2prot = {}
    for row in search_db(cursor,qry,False):
        [gene_id,  protein] =  row
        store_if(metacyc_gene2prot, gene_id, protein.replace(' ',''))
    return metacyc_gene2prot

#########################################
def mc_prot2gene_mapping(cursor):
    qry = 'select unique_id, gene  from metacyc_proteins where species="TAX-9606"'
    metacyc_prot2gene = {}
    for row in search_db(cursor,qry,False):
        [protein_id, gene] =  row
        if metacyc_prot2gene.has_key(protein_id):
            print "two genes for a protein? protein:", protein_id, "genes:", metacyc_prot2gene[protein_id], gene
            exit(1)
        else:
            metacyc_prot2gene[protein_id] = gene.replace(' ','')
    return metacyc_prot2gene

#########################################
def metacyc_uniprot2protein_mapping(metacyc_prot2gene, protein_id_translation, gene_id_translation):
    missed_keys = {}
    metacyc_uniprot2protein = {}
    #  manual corrections to metacyc
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

    return metacyc_uniprot2protein

##########################################
def get_mc_genes_from_mc_proteins(cursor):
    qry = "select distinct(gene) from metacyc_proteins where species='TAX-9606'"
    mc_genes = []
    for row in search_db(cursor, qry):
        for mc_gene in row[0].split(";"):
            mc_genes.append(mc_gene.replace(' ',''))
    return set(mc_genes)


##########################################
def create_metacyc_genes_entry(cursor, mc_gene, gene_id_translation):

    if not gene_id_translation.has_key(mc_gene) or len(gene_id_translation[mc_gene])<2:
        print "insufficient info - metacyc_genes entry not created"
        return
    [uniprot, common_name]  = gene_id_translation[mc_gene][:2]
    # what is the next available id?
    qry = "select max(id) from metacyc_genes"
    id = int(search_db(cursor, qry)[0][0])+1
    qry = "select unique_id from metacyc_proteins where gene='%s' " % mc_gene
    rows  = search_db(cursor, qry)
    if not rows:
        print "no product found in metacyc_proteins - metacyc_genes entry not created"
        return

    product = ";".join([row[0] for row in rows])
    qry = "insert into metacyc_genes "
    qry += "(id, common_name, unique_id, product) values "
    qry += "(%d, '%s', '%s', '%s') " %(id, common_name, mc_gene, product)
    search_db(cursor, qry, verbose=True)
    return
##########################################
def connect():
    development = gethostname()=='pegasus'
    if development:
        db = connect_to_mysql(user="cookiemonster", passwd=(os.environ['COOKIEMONSTER_PASSWORD']))
    else:
        db = connect_to_mysql(user="blimps", passwd=(os.environ['BLIMPS_DATABASE_PASSWORD']))
    if not db: exit(1)
    cursor = db.cursor()
    qry = 'set autocommit=1' # not sure why this has to be done explicitly - it should be the default
    search_db(cursor,qry,False)
    if development:
        switch_to_db(cursor, 'blimps_development')
    else:
        switch_to_db(cursor, 'blimps_production')
    return db, cursor

##########################################
def main():

    protein_id_translation = read_translation_table("/databases/meta_cyc/data/datfiles/protein-links.dat")
    gene_id_translation    = read_translation_table("/databases/meta_cyc/data/datfiles/gene-links.dat")
    # Weirdness: CKMT1A, CKMT1B  Two genes located near each other on chromosome 15 have been identified
    # which encode identical mitochondrial creatine kinase proteins, P12532
    # the genes: ENSG00000237289, ENSG00000223572
    # http://www.genecards.org/cgi-bin/carddisp.pl?gene=CKMT1A
    # similar thing for calmodulins 1, 2, 3
    # three genes, on completely different chromosomes (14, 2 an 19)
    # "encode an identical calcium binding protein which is one of the four subunits of phosphorylase kinase"
    # from http://www.genecards.org/cgi-bin/carddisp.pl?gene=CALM1

    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    db, cursor = connect()

    # direct mapping
    metacyc_prot2gene = mc_prot2gene_mapping(cursor)
    metacyc_gene2prot = mc_gene2prot_mapping(cursor)

    # inversion check:
    inversion_check(metacyc_gene2prot, metacyc_prot2gene, gene_id_translation)

    # roundabout mapping    # use uniprot ids to find metacyc_protein and metacyc_genes_ids

    metacyc_uniprot2protein = metacyc_uniprot2protein_mapping(metacyc_prot2gene, protein_id_translation, gene_id_translation)

    print " metacyc_prot2gene:", len(metacyc_prot2gene)
    print " metacyc_gene2prot:", len(metacyc_gene2prot)
    print " uniprot ids mapped to metacyc protein:", len(metacyc_uniprot2protein)
    print
    # scan through gene table
    # store them in gene table (add the columns to the table if needed)
    error_ct = 0
    qry = "select id, symbol, synonyms,  uniprot_ids from genes"
    #for row in search_db(cursor, qry):
    for row in []:
        [id, blimps_symbol, blimps_synonyms,  uniprot_ids] = row
        if not uniprot_ids: continue
        if not blimps_synonyms: blimps_synonyms = ''
        qry = "select id, unique_id, common_name from metacyc_genes where common_name= '%s'" % blimps_symbol
        rows = search_db(cursor, qry)
        if rows and len(rows)==1:
            qry = "update genes set metacyc_gene_id='%s' where id=%d " % (rows[0][1], id)
            search_db(cursor, qry, verbose=False)
            qry = "update metacyc_genes set gene_id='%s' where id=%d " % (id, rows[0][0])
            search_db(cursor, qry, verbose=False)

        elif uniprot_ids and uniprot_ids!="":
            # try searching by protein
            for uniprot_id in uniprot_ids.split ("|"):
                # metacyc proteins should be (more or less) metabolic enzymes only
                # the uniprot ids for the rest cannot be matched to metacyc
                if not metacyc_uniprot2protein.has_key(uniprot_id) : continue
                protein_ids = metacyc_uniprot2protein[uniprot_id].split(";")
                mc_gene_id = None
                for protein_id in protein_ids:
                    if not metacyc_prot2gene.has_key(protein_id):
                        print "no gene id for", protein_id
                        continue
                    if not mc_gene_id:
                        mc_gene_id =  metacyc_prot2gene[protein_id]
                        qry = "update genes set metacyc_gene_id='%s' where id=%d " % (mc_gene_id, id)
                        search_db(cursor, qry, verbose=False)
                        qry = "update metacyc_genes set gene_id=%d where unique_id='%s' " % (id, mc_gene_id)
                        search_db(cursor, qry, verbose=False)
                    elif mc_gene_id != metacyc_prot2gene[protein_id]:
                        print "gene id mismatch for", protein_ids

    # some mc_proteins are still unmapped
    if True: # this needs to be doen twice; the whole process needs to be cleaned up; this is a mess
        mc_genes = get_mc_genes_from_mc_proteins(cursor)

        qry = "select distinct(metacyc_gene_id) from genes where metacyc_gene_id is not null"
        blimps_genes = set([row[0] for row in search_db(cursor, qry)])
        print "human genes referenced in metacyc proteins:", len(mc_genes)
        print "blimps genes mapped to metacyc genes:", len(blimps_genes)
        print "mc genes not mapped:", len(mc_genes-blimps_genes)
        for mc_gene in mc_genes-blimps_genes:
            print "=== "  + mc_gene
        for mc_gene in mc_genes - blimps_genes:
            print "\n +++++++++++++ " + mc_gene
            qry = "select common_name from metacyc_genes where unique_id='%s'" % mc_gene
            ret = search_db(cursor, qry)
            if ret:
                if len(ret)==1:
                    symbol = ret[0][0]
                    if not symbol or symbol == "":
                        print mc_gene, " has empty string for common_name/symbol"
                        if gene_id_translation.has_key(mc_gene):
                            print "gene id translation: ", gene_id_translation[mc_gene]
                        else:
                            print "no gene id translation"
                    else:
                        qry = "select count(1) from genes where symbol='%s'" % symbol
                        count = int(search_db(cursor, qry, verbose=True)[0][0])
                        if count == 0:
                            print "no genes found with symbol", symbol
                            qry = "select count(1) from genes where synonyms like '%%%s%%'" % symbol
                            count = int(search_db(cursor, qry, verbose=True)[0][0])
                            print "no genes with", symbol, "as synonym: ", count
                            if count==1:
                                qry  = "update genes set metacyc_gene_id='%s' " % mc_gene
                                qry += "where synonyms like '%%%s%%'" % symbol
                                search_db(cursor, qry, verbose=True)
                        else:
                            qry = "update  genes set metacyc_gene_id='%s' where symbol='%s'" % (mc_gene, symbol)
                            search_db(cursor, qry, verbose=True)

            else:
                print qry
                if not ret:
                    print "no return for common name"
                    if gene_id_translation.has_key(mc_gene):
                        print "gene id translation: ", gene_id_translation[mc_gene]
                        # does this entry exist at all?
                        # example of nonexistent: HS04095, 'P21912', 'SDHB'
                        qry = "select count(1) from metacyc_genes where unique_id='%s'" % mc_gene
                        count = int(search_db(cursor, qry, verbose=True)[0][0])
                        if count==0: # not sure anymore how this happened
                            print "the entry", mc_gene, "does not exist ... creating"
                            create_metacyc_genes_entry(cursor, mc_gene, gene_id_translation)


                    else:
                        print " no gene id translation"
                else:
                    print "the return length", len(ret)

    cursor.close()
    db.close()
    return True


#########################################
if __name__ == '__main__':
    main()




