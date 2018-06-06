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
	target_vs_pct_idtty = find_pdb_ids_of_similar_seqs(sequence, 40,
									'./', qry_name, verbose=False)
	unique_pdbs = []
	for tgt, pct in target_vs_pct_idtty.iteritems():
		pdb_id = tgt[:4]
		if not pdb_id in unique_pdbs: unique_pdbs.append(pdb_id)

	for pdb_id in unique_pdbs:
		titles = pdb_titles(pdb_id)
		if len(titles)>1:
			print "multiple titles (?)"
			print "\n".join(titles)
			exit(1)
	    # not a good idea - some of the largest pieces of structure
		# were cy\srytallized with inhibitor
		#not_interesting = False
		#for smallmol in ['compound','peptid', 'inhibitor',
		#                 'nutlin', 'antagonist']:
		#	if smallmol in titles[0].lower(): not_interesting=True
		#if not_interesting: continue
		print "="*30
		print pdb_id, titles[0]
		for desc in mol_descriptions_from_pdb(pdb_id):
			print "\t", desc
			for chain in desc[1].split(","):
				chain_id = "_".join([pdb_id,chain])
				if chain_id in target_vs_pct_idtty.keys():
					print "\t\t", chain, target_vs_pct_idtty[chain_id]

#########################################
if __name__ == '__main__':
	main()
