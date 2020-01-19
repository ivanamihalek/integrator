#!/usr/bin/python

from integrator_utils.mysql import *

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
def single_nt_expansion(list_of_alts):

	all_divisible_by_3 = True
	for x in list_of_alts:
		if len(x)<2: continue
		x = x[1:]
		if x != x[0]*len(x): return False #this is not single nt
		if len(x)%3!=0: all_divisible_by_3 = False

	if all_divisible_by_3: return False # maybe this is 3nt expansion in disguise

	return True

##########################################
def main():

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
		# I am looking for triplets
		if len([x for x in list_of_alts if three_nt_change(len(x),ref_len)])==0: continue
		# I am looking for self similar triplets (with period of three, at least once
		if (not self_similar(ref,3)) and len([x for x in list_of_alts if self_similar(x,3)])==0: continue
		# I don't want single nucleotide expansions either
		if single_nt_expansion(list_of_alts): continue
		# try in the reverse direction too
		if single_nt_expansion([x[::-1] for x in list_of_alts]): continue
		# store all such cases, and see if some overlap and need to be merged
		max_reach = pos +max([ref_len]+ [len(x) for x in list_of_alts] ) - 1
		candidates.append([pos, ref, alt, var_counts, total_count, max_reach]) # I want the positions to remain sorted

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

	for cluster in clusters:
		if len(cluster)<2: continue
		print "********************"
		for candidate in cluster:
			print candidate
		print

	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914