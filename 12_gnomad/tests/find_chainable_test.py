#!/usr/bin/python

from integrator_utils.periodic_mods import *

#######################################
def find_comparable_patterns(raw_seqs, motifs):
	clusters = []
	raw_seqs.sort(lambda x, y: cmp(-len(x), -len(y)))
	for new_seq in raw_seqs:
		new_pattern = decomposition(new_seq, motifs)
		#print new_seq, new_pattern
		cluster_found = False
		for cluster in clusters:
			#patterns for all members of a cluster should be the same
			pattern = cluster[0]
			if equivalent(pattern,new_pattern):
				cluster.append(new_pattern)
				cluster_found = True
				break
		# if not seq placed, start new cluster
		if not cluster_found:
			clusters.append([new_pattern])

	return clusters

#########################################
def main():
	#raw_seqs = [u'AGAAGATGAT', u'AGAT', u'A', u'AGAAGAT', u'AGAA', u'AGATGATGAT', u'AGAAGATGATGAT', u'AGAAGATGATGATGAT', u'AGAAGAAT', u'AGAAGATGA']
	#motifs = ['GAT']
	#raw_seqs = [u'CCAGTCTTT', u'CGTCTTT', u'CCA', u'CCAGTCT']
	#motifs = ['T']
	#raw_seqs=['GCACCC','GCA']
    #motifs=['C']
	raw_seqs =  ['TGT','TT','TG','TGTT','TTT']
	motifs=['T']
	'''
	raw_seqs = ["CTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGCTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGC", "CTGCTGTTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGC",
				"CTGCTGTTGCTGCTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGTTGCTGC",
				"TGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGC",
				"CTGCTGTTGCTGTTGCTGT",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGTTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGTTGTTGCTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGTTGTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGCTGC",
				"CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGC"]

	motifs = ["TGC", "TGCTGT"]
	'''
	tot_cluster_length = 0
	for cluster in find_comparable_patterns(raw_seqs, motifs):
		tot_cluster_length += len(cluster)
		if len(cluster)<2: continue
		print "cluster:             "
		for pattern in cluster:
			print pattern
		for pattern in cluster:
			print to_string(pattern)
		print"--------------"
	print
	print len(raw_seqs), tot_cluster_length
	return

#########################################
if __name__ == '__main__':
	main()

