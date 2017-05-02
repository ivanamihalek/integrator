#!/usr/bin/python
from collections import deque

#######################################
def decompositions(string, motif):
	# decomposition is left_pad, multiplier (how many times does the motif repeat), right_pad
	l = len(motif)
	match_positions = deque()
	for i in range(len(string)):
		if string[i:i+l]==motif:
			match_positions.append(i)

        # the trivial decomposition, if there is no match inside
        if len(match_positions)==0:
	        decomps = [[string, 1, ""]]
                return decomps

        decomps = []
	multiplier = 0
	start      = -1
	while match_positions:
		mp = match_positions.popleft()
                if multiplier==0: start = mp
                multiplier += 1
		if len(match_positions)==0 or mp+l != match_positions[0]:
                   decomps.append([string[:start], multiplier, string[start+multiplier*l:] ])
                   multiplier = 0
                
	return decomps

#######################################
def find_chainable(raw_seqs, motifs):

	for motif in motifs:
		clusters = []
		cluster_pads = []
		# sort descending, by length
                raw_seqs.sort(lambda x, y: cmp(-len(x), -len(y)))
		for new_seq in raw_seqs:
			for decomposition in  decompositions(new_seq, motif):
 				[padding_left, multiplier, padding_right] = decomposition
				cluster_found = False
				for cluster in clusters:
					cluster_index = clusters.index(cluster)
                                        print  cluster_pads[cluster_index], [padding_left, padding_right]
					if cluster_pads[cluster_index] == [padding_left, padding_right]:
						cluster.append(new_seq)
						cluster_found = True
						break
				# if not seq placed, start new cluster
				if not cluster_found:
  					clusters.append([new_seq])
					cluster_pads.append([padding_left, padding_right])
		print "-------" , motif
		for cluster in clusters:
			if len(cluster)<2: continue
			# new_seq to cluster
			print "cluster:"
			cluster.sort(lambda x, y: cmp(len(x), len(y)))
			for seq in cluster: print seq
		print"--------------"
	print
	return

#########################################
def main():
    #raw_seqs = [u'AGAAGATGAT', u'AGAT', u'A', u'AGAAGAT', u'AGAA', u'AGATGATGAT', u'AGAAGATGATGAT', u'AGAAGATGATGATGAT', u'AGAAGAAT', u'AGAAGATGA']
    #motifs = ['GAT']
    raw_seqs = [u'CCAGTCTTT', u'CGTCTTT', u'CCA', u'CCAGTCT']
    motifs = ['T']
    find_chainable(raw_seqs, motifs)
    return

#########################################
if __name__ == '__main__':
	main()

