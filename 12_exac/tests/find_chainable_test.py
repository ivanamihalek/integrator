#!/usr/bin/python
from collections import deque


#######################################
class FiniteStateMachine:

	def __init__(self, string, motifs):
		self.string = string
		self.motifs = motifs
		self.motifs.sort(lambda x, y: cmp(-len(x), -len(y)))
		self.pattern = []
		self.no_motif_string = ""

	def find_pattern(self):
		self.no_match_state()
		return self.pattern

	def match_state(self, motif):
		if not motif or len(motif)==0: return
		count = 0
		while len(self.string) and self.string[:len(motif)]==motif:
			count += 1
			self.string = self.string[len(motif):]
		if count > 0:  self.pattern.append([motif,count])
		self.no_match_state()
		return

	def no_match_state (self):
		matched_motif = None
		while len(self.string) and not matched_motif:
			mms = filter (lambda mm: self.string[:len(mm)]==mm,  self.motifs)
			if len(mms)==0: # no motif
				self.no_motif_string +=  self.string[0]
				self.string = self.string[1:]
			else:
				matched_motif = mms[0]

		self.emit_no_motif()
		self.match_state(matched_motif)
		return

	def emit_no_motif(self):
		if len(self.no_motif_string) > 0: self.pattern.append([self.no_motif_string, 1])
		self.no_motif_string = ""
		return

#######################################
def decomposition(string, motifs):
	fsm = FiniteStateMachine(string, motifs)
	return fsm.find_pattern()


#######################################
def equivalent(p1, p2):
	# pattern is a list of pairs [motif, number_of_repeats]
	# two patterns are equivalent if the order of motifs is the same,
	# irrespective of the number of times they repeat
	number_of_motifs = len(p1)
	if number_of_motifs != len(p2): return False
	matching_motifs = filter(lambda i: p1[i][0] == p2[i][0], range(number_of_motifs))
	return len(matching_motifs) == number_of_motifs

def to_string(pattern):
	string = ""
	for [motif,number_of_repeats] in pattern:
		string += motif*number_of_repeats
	return string

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
    #raw_seqs =  ['TGT','TT','TG','TGTT','TTT']
    #motifs=['T']

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

