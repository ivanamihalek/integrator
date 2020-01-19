#!/usr/bin/python

from mysql import *

def three_nt_change(len_a,len_b):
	if len_a==len_b: return False
	if (len_a-len_b)%3==0: return True
	return False

def self_similar (input_string, period):
	s1 = input_string[1:]
	if len(s1)<=period: return False
	# circular shift
	s2 = s1[period:] + s1[:period]
	common_ntuples = 0
	#for i in range(0, len(s1)-period, period):
	# for now, look only for immediate repetition
	i = 0
	if s1[i:i+period] == s2[i:i+period]: common_ntuples += 1

	return common_ntuples>0

##########################################
def minimal_motif(ym):
	for i in range(1, len(ym) / 2 + 1):
		circular = ym[i:] + ym[:i]
		if circular == ym:
			return ym[:i]  # this is periodic with period i
	return ym

###############
def find_motif_in_pair(alt, ref,  find_minimal = True):
	if len(alt) > len(ref):
		[A, B] = [alt, ref]
	else:
		[A, B] = [ref, alt]
	if A[:len(B)] != B: return None  # this is not what I'm looking for
	ym = A[len(B):]
	# the simplest  possibility: x==y
	if len(B) >= len(ym) and B[-len(ym):] == ym:
		if not find_minimal or len(ym) == 1:
			return ym
		else:
			return minimal_motif(ym)
	mm = minimal_motif(ym)
	if len(mm) < len(ym):
		return mm

	return None



###############
def find_motif_in_variant(v, find_minimal = True):

	# for each pair reference and alt seq  - take the A to be longer of the two, and B the shorter
	# assume A = a + xm +ym, and B = a + xm
	# where a is an aperiodic piece, and m is a periodically repeating motif
	# (repeating x times in A and y times in B)
	# ym  = A - B
	# the period y to be determined by circular shifting
	[pos, ref, alts, var_counts, total_count, max_reach] = v
	for alt in alts.split(","):
		mm = find_motif_in_pair(alt, ref, find_minimal)
		if not mm:
			continue
		else:
			return mm
	return None

##########################################
def has_repeating_motif(v):
	if find_motif_in_variant(v, find_minimal=False):
		return True
	return False

##########################################
def periodicity_found(cluster):
	for variant in cluster:
		if has_repeating_motif(variant):return True
	# no variant has a repeating motif
	return False

##########################################
def main():

	if False:
		v = [16269823L, 'A', 'T,AAAAAAT', '13,3', 69068L, 16269829L]
		find_motif_in_variant(v, find_minimal=False)
		find_motif_in_variant(v, find_minimal=True)
		exit()

	db, cursor = connect()
	chrom = "22"
	table = "exac_freqs_chr_" + chrom
	print "*"*20
	print table
	print "number of variants:", search_db(cursor, "select count(1) from %s" % table)[0][0]
	qry  = "select count(1) from %s " % table
	qry += "where char_length(reference)=1 and char_length(variants)=1"
	print "simple SNPs",  search_db(cursor, qry)[0][0]
	print
	print "complex variants"
	qry  = "select * from %s " % table
	qry += "where char_length(reference)!=1 or char_length(variants)!=1"

	candidates = []
	for variant in search_db(cursor, qry):
		[pos, ref, alt, var_counts, total_count] = variant
		# from these further remove cases where the variant field is a list of SNPs
		list_of_alts = alt.split(",")
		ref_len = len(ref)
		if ref_len==1 and ("," in alt) and (len(alt)+1)==2*len(list_of_alts): continue
		max_reach = pos +max([ref_len]+ [len(x) for x in list_of_alts] ) - 1
		candidates.append([pos, ref, alt, var_counts, total_count, max_reach]) # I want the positions to remain sorted

	print "Done scanning. Looking for clusters."
	clusters = []
	for candidate in candidates:
		cluster_found = False
		for cluster in clusters:
			if len([x for x in cluster if x[0] <= candidate[0] <= x[-1]]):
				cluster.append(candidate)
				cluster_found = True
				break
		if not cluster_found: # start new cluster - cluster is a list of candidates
			clusters.append([candidate])

	print "Done clustering. Looking for repeats."

	# TODO: even without periodicity some clusters seem to describe related events

	clusters_w_periodicity = []
	for cluster in clusters:
		if len(cluster)<2: continue
		if not periodicity_found(cluster): continue
		clusters_w_periodicity.append(cluster)

	for cluster in clusters_w_periodicity:
		print "********************"
		smallest_pos = min ([x[0] for x in cluster])
		for v in cluster:
			[pos, ref, alts, var_counts, total_count, max_reach] = v
			motif = find_motif_in_variant(v)
			print " %12d   %s           %s " % (pos, " "*(pos-smallest_pos) + ref, alts),
			if motif: print "  motif : %s" % motif,
			print
		print

	print len(clusters)
	print
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

