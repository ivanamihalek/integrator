#!/usr/bin/python
from collections import deque

#######################################
def decompositions(string, motif):

        # the trivial decomposition, if there is no match inside
        decomps = [[string, 0, ""]]
        
	# decomposition is left_pad, multiplier (how many times does the motif repeat), right_pad
	l = len(motif)
	match_positions = deque()
	for i in range(len(string)):
		if string[i:i+l]==motif:
			match_positions.append(i)

        if len(match_positions)==0:
               return decomps

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
def equivalent_padding(padding1, padding2):

        #return padding1==padding2

        [left1, right1] = padding1
        [left2, right2] = padding2

        
        left_ok = left1==left2

        if right1==right2:
                right_ok = True
        else:
                shorter = len(right1)
                if len(right1)>len(right2):  shorter = len(right2)

                right_ok = (shorter>2) and right1[:shorter]== right2[:shorter]
       
        return left_ok and right_ok

#######################################
def find_chainable(raw_seqs, motifs):
	clusters = {}
	cluster_pads = {}
	for motif in motifs:
		clusters[motif] = []
		cluster_pads[motif] = []
                raw_seqs.sort(lambda x, y: cmp(-len(x), -len(y)))
		for new_seq in raw_seqs:
			# special case: motif is an insert in the new seq:
			if not motif in new_seq:
				for cluster_index in range(len(clusters[motif])):
					cluster = clusters[motif][cluster_index]
					if "".join(cluster_pads[motif][cluster_index]) == new_seq:
						cluster.append(0)
			else:
                                for decomposition in  decompositions(new_seq, motif):
                                        [padding_left, multiplier, padding_right] = decomposition
                                        print new_seq, decomposition
                                        cluster_found = False
                                        for cluster_index in range(len(clusters[motif])):

                                                cluster = clusters[motif][cluster_index]
                                                #print  cluster_pads[cluster_index], [padding_left, padding_right]
                                                if equivalent_padding (cluster_pads[motif][cluster_index], [padding_left, padding_right]):
                                                        cluster.append(multiplier)
                                                        cluster_found = True
                                                        break
                                        # if not seq placed, start new cluster
                                        if not cluster_found:
                                                clusters[motif].append([multiplier])
                                                cluster_pads[motif].append([padding_left, padding_right])
                
	return clusters, cluster_pads

#########################################
def main():
    #raw_seqs = [u'AGAAGATGAT', u'AGAT', u'A', u'AGAAGAT', u'AGAA', u'AGATGATGAT', u'AGAAGATGATGAT', u'AGAAGATGATGATGAT', u'AGAAGAAT', u'AGAAGATGA']
    #motifs = ['GAT']
    raw_seqs = [u'CCAGTCTTT', u'CGTCTTT', u'CCA', u'CCAGTCT']
    motifs = ['T']
    #raw_seqs=['GCACCC','GCA']
    #motifs=['C']
    #raw_seqs =  ['TGT','TT','TG','TGTT','TTT']
    #motifs=['T']
    clusters_with_motif, cluster_pads = find_chainable(raw_seqs, motifs)
    if True:
	for motif, clusters in clusters_with_motif.iteritems():
		print "------- motif: ", motif
		for cluster_index in range(len(clusters)):
			cluster = clusters[cluster_index]

			if len(cluster)<2: continue
			print "cluster:             ", cluster_pads[motif][cluster_index]
			cluster.sort()
			left_pad  = cluster_pads[motif][cluster_index][0]
			right_pad = cluster_pads[motif][cluster_index][1]
			for multiplier in cluster:
				raw_seq = left_pad + motif*multiplier + right_pad
				print multiplier, raw_seq
		print"--------------"
	print


    
    return

#########################################
if __name__ == '__main__':
	main()

