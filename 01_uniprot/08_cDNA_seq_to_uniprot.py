#! /usr/bin/python

# getting the cdna
# ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cds

from integrator_utils.mysql import *
import os,  subprocess
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
###################################################
def cleanseq(seq):
	clean = ""
	for line in seq.split("\n"):
		line.replace(" ","")
		if len(line)==0: continue
		if line[0]=='>': continue
		line.replace(" ","")
		clean += line
	return clean
###################################################
def find_relative_transl_start(hg19_exon_starts, hg19_exon_ends, hg19_cds_start):
	rel = 0
	starts = [int(x) for x in hg19_exon_starts.split(",")]
	ends   = [int(x) for x in hg19_exon_ends.split(",")]
	if (len(starts)!=len(ends)):
		print "start end elngth mimsmatch"
		exit(1)

	for i in range(len(starts)):
		if starts[i] <= hg19_cds_start < ends[i]:
			rel += hg19_cds_start - starts[i]
			break
		else:
			rel += ends[i] - starts[i]
	return rel

###################################################
def get_seq_from_cds(blastextract, cdsdb, ensid):
	cmd = "{} -db {} -entry {}".format(blastextract, cdsdb,ensid)
	cdnaseq = subprocess.check_output(cmd, shell=True)
	if 'Error' in cdnaseq: return None
	cleandna = cleanseq(cdnaseq)
	new_biopython_seq = Seq(cleandna, unambiguous_dna)
	translated_seq = new_biopython_seq.translate().split("*")[0]
	return cleandna, translated_seq


###################################################
def get_seq_from_cdna(blastextract, cdnadb, ensid, hg19_exon_starts, hg19_exon_ends, hg19_cds_start, aa_sequence):
	cmd = "{} -db {} -entry {}".format(blastextract, cdnadb,ensid)
	cdnaseq = subprocess.check_output(cmd, shell=True)
	if 'Error' in cdnaseq: return None
	cleandna = cleanseq(cdnaseq)
	transl_start = find_relative_transl_start(hg19_exon_starts, hg19_exon_ends, hg19_cds_start)
	new_biopython_seq = Seq(cleandna[transl_start+1:], unambiguous_dna)
	seq = new_biopython_seq.translate().split("*")[0]
	if seq==aa_sequence:
		print "cdna translation ok"
		return seq
	# failure - try all reading frames
	for offset in [0,1,2]:
		new_biopython_seq = Seq(cleandna[offset:], unambiguous_dna)
		if aa_sequence in new_biopython_seq:
			print "found in "
			print new_biopython_seq
			exit()
	return None

###################################################
def find_coding_dna_sequence (cursor, blastextract,  cdsdb,cdnadb, ensembl_gene_id):
	switch_to_db(cursor, 'blimps_development')
	qry  = "select uniprot_id from uniprot_basic_infos "
	qry += "where ensembl_gene_id like '%%%s%%'" % ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret:
		print ensembl_gene_id, "not found in uniprot_basic_infos"
		exit()
	switch_to_db(cursor, 'monogenic_development')
	for line in ret:
		uniprot_id = line[0]
		qry  = "select ensembl_transcript_id, sequence from uniprot_seqs where uniprot_id='%s'" % uniprot_id
		ret2 =  search_db(cursor,qry)
		if not ret2:
			print uniprot_id, ensembl_gene_id, "not found"
			continue
		for line2 in ret2:
			[ensembl_transcript_id, uniprot_sequence] = line2
			# traascript id has soem sub-something indicated by .digit after the transcript number
			cmd = "grep %s %s | awk '{print $1}' | sed 's/>//g'" % (ensembl_transcript_id, cdsdb)
			grepret = subprocess.check_output(cmd, shell=True).rstrip().split("\n")
			for ensid in grepret:
				if len(ensid) == 0: continue
				cds, cds_translation = get_seq_from_cds(blastextract, cdsdb, ensid)
				if not cds_translation:
					print "no cds translation for ", ensid
				qry  = "update uniprot_seqs set ensembl_coding_sequence='%s'" % cds
				qry += "where uniprot_id= '%s' " % uniprot_id
				search_db(cursor,qry, verbose=False)
				if  uniprot_sequence != cds_translation:
					qry  = "update uniprot_seqs set ensembl_cds_translation='%s'" % cds_translation
					qry += "where uniprot_id= '%s' " % uniprot_id
					search_db(cursor,qry, verbose=False)
	return

###################################################
def main():

	blastextract = "/usr/local/bin/blastdbcmd"
	cdsdb = "/databases/ensembl/hg19/cds/Homo_sapiens.GRCh38.cds.all.fa"
	cdnadb = "/databases/ensembl/hg19/cdnaseq/Homo_sapiens.GRCh38.cdna.all.fa"
	for dependency in [blastextract, cdsdb, cdnadb]:
		if not os.path.exists(dependency):
			print dependency, "not found"
			exit(1)

	db, cursor = connect()
	# make sure that the uniprot_seqs table has cds columns
	for col in  ['ensembl_coding_sequence', 'ensembl_cds_translation']:
		if 'coding' in col:
			coltype= "mediumtext"
		else:
			coltype= "text"
		check_or_make_column (cursor, 'monogenic_development', 'uniprot_seqs', col, coltype)

	# limit ourselves to iems related genes
	qry = "select ensembl_gene_id from omim_genemaps where inborn_error_of_metabolism=1"
	ret = search_db(cursor, qry)
	for line in ret:
		ensembl_gene_id = line[0]
		find_coding_dna_sequence (cursor, blastextract,  cdsdb,cdnadb, ensembl_gene_id)



	return


#########################################
if __name__ == '__main__':
	main()
