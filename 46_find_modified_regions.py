#!/usr/bin/python

from integrator_utils.mysql import *

#########################################
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
		max_reach = pos +max([ref_len]+ [len(x) for x in list_of_alts] ) - 1
		candidates.append([pos, ref, alt, var_counts, total_count, max_reach]) # I want the positions to remain sorted

	print "Done scanning. Looking for clusters."
	clusters = []
	for candidate in candidates:
		cluster_found = False
		for cluster in clusters:
			if len([x for x in cluster if x[0] <= candidate[0] <= max(x[0]+3,x[-1]) ]):
				cluster.append(candidate)
				cluster_found = True
				break
		if not cluster_found: # start new cluster - cluster is a list of candidates
			clusters.append([candidate])

	print "Done clustering. Looking for repeats."

	ct = 0
	for cluster in clusters:
		if len(cluster)<2: continue
		print "********************"
		ct += 1
		smallest_pos = min ([ x[0] for x in cluster])
		largest_pos  = max ([x[-1] for x in cluster])
		for v in cluster:
			[pos, ref, alts, var_counts, total_count, max_reach] = v
			for sequence in [ref] + alts.split(","):
				print " %12d   %s" % (pos, "-"*(pos-smallest_pos) + sequence + "-"*(largest_pos - (pos+len(sequence)-1)))
		print

	print ct
	print
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914