#!/usr/bin/python
#
#
#

from integrator_utils.python.mysql     import *
from integrator_utils.python.restfuls  import *
from integrator_utils.python.pdb_tools import *
import integrator_utils.python.pdb_constants as pdbc

#########################################
def readseq(seqfile):
	sequence = ""
	seqfile_handle = open(sys.argv[1],"r")
	for line in seqfile_handle:
		if line[0]!='>':  sequence += line.rstrip()
	seqfile_handle.close()
	return sequence
#########################################
def main():

	################# assorted input checks #######################
	if not len(sys.argv)>1:
		print "Usage: %s <seq_file (fasta)> " % sys.argv[0]
		exit(0)
	pdb_check_resources()
	seqfile_name = sys.argv[1]
	if seqfile_name[-6:]!='.fasta':
		print 'Expected .fasta extension in the sequence file name'
		exit(1)
	qry_name = seqfile_name[:-6]


	sequence = readseq(seqfile_name)
	find_pdb_ids_of_similar_seqs(sequence, 40, './', qry_name, verbose=True)


#########################################
if __name__ == '__main__':
	main()
