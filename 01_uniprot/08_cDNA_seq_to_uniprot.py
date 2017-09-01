#! /usr/bin/python

# getting the cdna
# ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cds

from integrator_utils.mysql import *
import os,  subprocess, requests, re
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
def get_coding_ranges(hg19_exon_starts, hg19_exon_ends, hg19_cds_start, hg19_cds_end):
	# All coordinates in these tables are half-open zero-based -- note though that the BROWSER AND EXAC ARE NOT 0-BASED!
	# This means that the first 100 bases of a chromosome are represented as [0,100), i.e. 0-99.
	# The second 100 bases are represented as [100,200), i.e. 100-199, and so forth.
	# An advantage of half-open coordinate ranges is that
	# the length can be obtained by simply subtracting the start from the end.

	exon_starts = [int(i) for i in hg19_exon_starts.split(",")]
	exon_ends   = [int(i) for i in hg19_exon_ends.split(",")]
	cds_start   = int(hg19_cds_start)
	cds_end     = int(hg19_cds_end)
	coding_exons = []
	tot_number_of_exons = len(exon_starts)
	for i in range(tot_number_of_exons):
		s = exon_starts[i]
		e = exon_ends[i]
		if cds_start>e: continue
		coding_start = s
		if cds_start>=s:
			coding_start = cds_start
		if cds_end<s: break
		coding_end = e
		if cds_end<=e:
			coding_end = cds_end
		coding_exons.append([coding_start+1,coding_end+1])

	return coding_exons

#########################################
def get_region_from_das(assembly, chrom, start, end):
	das_request  = "http://genome.ucsc.edu/cgi-bin/das/%s/" % assembly
	das_request += "dna?segment=chr%s:%s,%s" %(chrom, start, end)

	ret = requests.get(das_request).text.lower().replace("\n","")
	pattern = re.compile("<dna.*?>([\w\s]*)</dna>")
	matches = re.findall(pattern, ret)
	dna = ""
	for m in matches: dna += m.replace(" ","")
	return dna


###################################################
def sequence_from_das(cursor, uniprot_id):
	# Example: ENST00000377332, the ensembl translation doe not correspond to the one
	# given in ucsc (as enembl seqeunce)
	qry = "select * from monogenic_development.uniprot_seqs where uniprot_id= '%s' " % uniprot_id

	[uniprot_id, ensembl_transcript_id, ensembl_protein_id, chrom, strand,\
	hg19_exon_starts, hg19_exon_ends, hg19_cds_start, hg19_cds_end, \
	sequence, cons_in_verts_map_to_uni, ensembl_coding_sequence, \
	ensembl_cds_translation, uni2ens_cigar, cons_in_verts_map_to_ens] = search_db(cursor,qry)[0]

	coding_exons = get_coding_ranges(hg19_exon_starts, hg19_exon_ends, hg19_cds_start, hg19_cds_end)
	reconstructed_cds = ""

	for coding_range in coding_exons:
		reconstructed_cds += get_region_from_das('hg19', chrom, coding_range[0], coding_range[1]-1)
	new_biopython_seq = Seq(reconstructed_cds, unambiguous_dna)
	if strand=="-": new_biopython_seq = new_biopython_seq.reverse_complement()

	translated_seq = new_biopython_seq.translate().split("*")[0]

	return str(new_biopython_seq), str(translated_seq)

###################################################
def find_coding_dna_sequence (cursor,ensembl_gene_id):
	switch_to_db(cursor, 'blimps_development')
	qry  = "select uniprot_id from uniprot_basic_infos "
	qry += "where ensembl_gene_id like '%%%s%%'" % ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret:
		print ensembl_gene_id, "not found in uniprot_basic_infos"
		exit()
	if len(ret)>1:
		print "more than one uniprot entry for ensembl ",  ensembl_gene_id
		print ret
		print "while possible in principle, not equipped to deal with it here"
		exit()
	uniprot_id = ret[0][0]
	print uniprot_id

	switch_to_db(cursor, 'monogenic_development')
	qry  = "select  sequence from uniprot_seqs where uniprot_id='%s'" % uniprot_id
	ret =  search_db(cursor,qry)
	if not ret:
		print uniprot_id, "sequence not found"
		exit()

	for line in ret:
		uniprot_sequence = line[0]
		cds, cds_translation = sequence_from_das(cursor, uniprot_id)
		if not cds_translation:
			print "no cds translation for ", uniprot_id
		qry  = "update uniprot_seqs set ensembl_coding_sequence='%s'" % cds
		qry += "where uniprot_id= '%s' " % uniprot_id
		search_db(cursor,qry, verbose=False)
		print uniprot_sequence
		print '--------------------'
		print cds_translation

		if  uniprot_sequence == cds_translation:
			print "uniprot_sequence can be reconstructed from ensembl "
			qry  = "update uniprot_seqs  set ensembl_cds_translation=null "
			qry += "where uniprot_id= '%s' " % uniprot_id
			search_db(cursor,qry, verbose=False)
			in_agreement = True
			break
		else:
			print "uniprot_sequence does not equal cds translation - storing ensembl translation "
			qry  = "update uniprot_seqs  set ensembl_cds_translation='%s' " % cds_translation
			qry += "where uniprot_id= '%s' " % uniprot_id
			search_db(cursor,qry, verbose=False)
			in_agreement = False

		# in either case, set the alignemnt field to null because it needs reworking
		qry  = "update uniprot_seqs  set uni2ens_cigar=null "
		qry += "where uniprot_id= '%s' " % uniprot_id
		search_db(cursor,qry, verbose=False)
		print

	return in_agreement

###################################################
def main():

	db, cursor = connect()
	# make sure that the uniprot_seqs table has cds columns
	for col in  ['ensembl_coding_sequence', 'ensembl_cds_translation']:
		if 'coding' in col:
			coltype= "mediumtext"
		else:
			coltype= "text"
		check_or_make_column (cursor, 'monogenic_development', 'uniprot_seqs', col, coltype)

	# limit ourselves to iems related genes
	qry = "select ensembl_gene_id from blimps_development.omim_genemaps where inborn_error_of_metabolism=1"
	ensids  =[line[0] for line in  search_db(cursor, qry)]
	ensids = ['ENSG00000178057']
	total = 0
	in_agreement = 0
	for  ensembl_gene_id in ensids:
		print
		print ensembl_gene_id
		uniprot_and_ensembl_agree = find_coding_dna_sequence (cursor,  ensembl_gene_id)
		total += 1
		if uniprot_and_ensembl_agree: in_agreement+= 1

	print
	print "total IEM genes:", total
	print "uniprot and ensembl seq in agreement:", in_agreement
	print


	return


#########################################
if __name__ == '__main__':
	main()
