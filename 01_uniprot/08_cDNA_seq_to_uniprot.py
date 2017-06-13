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
	seq = new_biopython_seq.translate().split("*")[0]
	return seq

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
		return seq
	# failure - try all reading frames

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
	switch_to_db(cursor, 'monogenic_development')
	qry = "select * from uniprot_seqs"
	ret = search_db(cursor,qry)
	for line in ret:
		[uniprot_id, ensembl_transcript_id, ensembl_protein_id, chrom, strand,
		 hg19_exon_starts, hg19_exon_ends, hg19_cds_start, hg19_cds_end, sequence, conservation] = line
		# traascript id has soem sub-something indicated by .digit after the transcript number
		cmd = "grep %s %s | awk '{print $1}' | sed 's/>//g'" % (ensembl_transcript_id, cdsdb)
		grepret = subprocess.check_output(cmd, shell=True).rstrip().split("\n")
		for ensid in grepret:
			if len(ensid) == 0: continue
			cds_translation = get_seq_from_cds(blastextract, cdsdb, ensid)
			if  sequence != cds_translation:
				print
				print uniprot_id, ensembl_transcript_id
				print cds_translation
				print "---------------------------------"
				print sequence
				print grepret
				from_cdna = get_seq_from_cdna(blastextract, cdnadb, ensid, hg19_exon_starts, hg19_exon_ends, hg19_cds_start)
				print "---------------------------------"
				print from_cdna
				exit()

	return


#########################################
if __name__ == '__main__':
	main()
