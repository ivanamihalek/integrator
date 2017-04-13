#!/usr/bin/python

from integrator_utils.mysql import *

#########################################
def insert (alignment, s, relative_position, ins):

	# we are inserting in sequence s, at  relative_position
	for i in range(len(ins)): # position in the insert
		seq_pos  = relative_position + i
		almt_pos  = -1
		nongap_ct = 0
		while  almt_pos< len(alignment[s])-1 and nongap_ct < seq_pos:
			almt_pos += 1
			if alignment[s][almt_pos] != "-": nongap_ct += 1
		# insert at the following position:
		almt_pos += 1
		if almt_pos>=len(alignment[s]) or  alignment[s][almt_pos] != "-":
			for s2 in range(len(alignment)):
				new_seq = alignment[s2][:almt_pos]
				if s2==s:
					new_seq += ins[i]
				else:
					new_seq += "-"
				new_seq += alignment[s2][almt_pos:]
				alignment[s2] = new_seq
		relative_position  += 1
	return

#########################################
def reconstruct_alignment(cluster):
	smallest_ref_pos = min ([x[0] for x in cluster])
	upper_bound_ref = max ([x[0] + len(x[1]) for x in cluster])
	refseq = "-" * (upper_bound_ref - smallest_ref_pos)
	for v in cluster:
		[pos, ref, alts, var_counts, total_count, max_reach] = v
		relative_pos = pos - smallest_ref_pos
		new_seq = refseq[:relative_pos] + ref + refseq[relative_pos + len(ref):]
		refseq = new_seq
	alignment = []
	seqinfo   = []
	# formal "variant" to hold reference sequence
	alt = [smallest_ref_pos, refseq[0], refseq[0], 1, 1, 1]
	alignment.append(refseq)
	seqinfo.append(alt)
	for v in cluster:
		[pos, ref, alts, var_counts, total_count, max_reach] = v
		relative_pos = pos - smallest_ref_pos
		counts = var_counts.split(",")
		alternatives = alts.split(",")
		a = -1 #  index of the alternative sequence fo this position
		for sequence in alternatives:
			a += 1
			alt = [pos, ref, sequence, int(counts[a]), total_count, len(ref)]
			seqinfo.append(alt)
			refseq = alignment[0]
			if len(sequence) <= len(ref):
				seq_padded   = sequence + "-"*(len(ref)-len(sequence)) # negative padding == no padding
				modified_seq = refseq[:relative_pos] + seq_padded + refseq[relative_pos + len(ref):]
				alignment.append(modified_seq)
			else:
				# append the reference sequence, than put an insert
				alignment.append(refseq)
				s = len(alignment)-1 # the index of our sequence
				#insert (alignment,  s, relative_pos+len(ref), sequence[len(ref):])


	for s in range(len(alignment)):
		modified_seq = alignment[s]
		[pos, ref, sequence, count, total_count, max_reach] = seqinfo[s]
		print " %12d  %3d  %s" % (pos, pos-smallest_ref_pos+1, modified_seq),
		print "       %7d  %7d    %s    %s" % (count, total_count, ref[:10], sequence[:10])

	return

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
	for candidate in candidates[0:len(candidates)/10]:
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
		reconstruct_alignment(cluster)

	print ct
	print
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914