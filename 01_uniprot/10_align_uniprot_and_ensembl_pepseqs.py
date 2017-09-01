#! /usr/bin/python
from integrator_utils.mysql import *
from Bio.Align.Applications import MuscleCommandline
###################################################
def main():

	db, cursor = connect()
	switch_to_db(cursor,"monogenic_development")

	check_or_make_column (cursor, 'monogenic_development', 'uniprot_seqs', 'uni2ens_cigar','text')

	qry ="select uniprot_id, sequence, ensembl_cds_translation from uniprot_seqs where ensembl_cds_translation !=''"
	for line in  search_db(cursor,qry):
		[uniprot_id, sequence, ensembl_cds_translation] = line
		outf = open("tmp.fasta","w")
		outf.write(">{}\n{}\n>{}\n{}\n".format('uni',sequence,'ens',ensembl_cds_translation))
		outf.close()
		afa_file = "{}.afa".format(uniprot_id)
		cline = MuscleCommandline(input='tmp.fasta', out=afa_file,  verbose=True)
		stdout, stderr = cline()
		inf = open(afa_file,"r")
		aln_seq = {}
		for name in ['uni','ens']: aln_seq[name]=""
		for line in inf:
			line = line.rstrip()
			if line[0]==">":
				name = line[1:]
			else:
				aln_seq[name] += line
		inf.close()
		old_operation = None
		operation_ct = 0
		cigar = ""
		for i in range(len(aln_seq['uni'])):
			if aln_seq['uni'][i]==aln_seq['ens'][i]:
				operation = 'M'
			elif aln_seq['uni'][i]=='-':
				operation = 'D'
			elif aln_seq['ens'][i]=='-':
				operation = 'I'
			else:
				operation = 'X'
			if old_operation and old_operation!=operation:
				cigar += "{}{}".format(old_operation, operation_ct)
				operation_ct = 0
			old_operation = operation
			operation_ct += 1
		cigar += "{}{}".format(operation, operation_ct)
		print uniprot_id,cigar
		qry =  "update uniprot_seqs set uni2ens_cigar='%s' " % cigar
		qry += "where uniprot_id='%s'" % uniprot_id
		search_db(cursor,qry, verbose=True)

	cursor.close()
	db.close()
	return


#########################################
if __name__ == '__main__':
	main()
