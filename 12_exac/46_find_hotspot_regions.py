#!/usr/bin/python

# getting the  genomic region from UCSC DAS server
# https://github.com/ivanamihalek/progesterone/blob/master/PR_binding_sites/scripts/ruby_modules/httputils.rb
# python http request http://docs.python-requests.org/en/master/

from integrator_utils.mysql import *
from integrator_utils.periodic_mods import *

import requests # http  requests
import re
from collections import deque

assembly = "hg19"

#########################################
def to_almt_pos(sequence, nongapped_pos):
	almt_pos  = -1
	nongap_ct = -1
	for character in sequence:
		almt_pos += 1
		if character!="-": nongap_ct +=1
		if nongap_ct >= nongapped_pos: break

	if nongap_ct<nongapped_pos: # appending to the end
		almt_pos = len(sequence)
	return almt_pos

#########################################
def insert (alignment, s, reference_position, ins):
	# insert ins at  reference_pos
	for s2 in range(len(alignment)):
		aux = alignment[s2]
		if s2 == s:
			alignment[s2] = aux[:reference_position] + ins + aux[reference_position:]
		else:
			alignment[s2] = aux[:reference_position] + "-"*len(ins) + aux[reference_position:]
	return
#########################################
def replace_piece (refseq, relative_pos, newseq):
	modified_seq = refseq
	for sp in range(len(newseq)):
		almt_pos = to_almt_pos(refseq, relative_pos+sp)
		modified_seq = modified_seq[:almt_pos] + newseq[sp] + modified_seq[almt_pos+1:]
	return modified_seq

#########################################
def add_modified_seq(alignment, sequence, ref, relative_pos):
	refseq = alignment[0]
	if len(sequence) <= len(ref):
		seq_padded = sequence + "-" * (len(ref) - len(sequence))  # negative padding == no padding
		modified_seq = replace_piece(refseq, relative_pos, seq_padded)
		alignment.append(modified_seq)
	else:
		# append the reference sequence, with the original seq replaced with the
		# noninserting part of the new sequence
		modified_seq = replace_piece(refseq, relative_pos, sequence[:len(ref)])
		alignment.append(modified_seq)
		s = len(alignment) - 1  # the index of our sequence in the alignment
		# now add the inserting part in our sequence and gaps everywhere else
		almt_pos = to_almt_pos(refseq, relative_pos + len(ref))
		insert(alignment, s, almt_pos, sequence[len(ref):])

##################################################
def	number_of_appearances(motif, string):
	l = len(motif)
	count = 0
	start_position = -1
	for i in range(len(string)-l+1):
		if string[i:i+l]==motif:
			count += 1
			start_position = i
			break
	if start_position<0:
		return 0
	for i in range(start_position+l,len(string)-l+1,l):
		if string[i:i+l]==motif:
			count += 1
		else: # it has to be consecutive
			break
	return count

#######################################
def snip(string):
	r = string
	if len(string) > 10: r = string[:10] + "..."
	return r

#######################################
def decompositions(string, motif):
	# decomposition is left_pad, multiplier (how many times does the motif repeat), right_pad
	l = len(motif)
	match_positions = deque()
	for i in range(len(string)):
		if string[i:i+l]==motif:
			match_positions.append(i)

	if len(match_positions)==0:
		# the trivial decomposition, if there is no match inside
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
	print raw_seqs
	print motifs
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
def reconstruct_alignment(chrom, cluster):


	periodic = periodicity_found(cluster)
	if not periodic: return

	print "\n********************"
	print "periodic"

	smallest_ref_pos = min ([x[0] for x in cluster])
	upper_bound_ref = max ([x[0] + len(x[1]) -1 for x in cluster])
	original_region = get_region_from_das(assembly, chrom, smallest_ref_pos, upper_bound_ref).upper()
	refseq = original_region
	for v in cluster:
		[pos, ref, alts, var_counts, total_count, max_reach] = v
		relative_pos = pos - smallest_ref_pos
		new_seq = refseq[:relative_pos] + ref + refseq[relative_pos + len(ref):]
		refseq = new_seq
	print " %12d  %3d  %s" %  (smallest_ref_pos, 1, refseq)
	print "     -------------------------------------------"
	alignment = []
	seqinfo   = []
	motifs    = []
	# formal "variant" to hold reference sequence
	alt = [smallest_ref_pos, refseq[0], refseq[0], 1, 1, 1, None]
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
			motif = None
			if periodic: motif = find_motif_in_pair(sequence, ref)
			if motif and not motif in motifs: motifs.append(motif)
			alt = [pos, ref, sequence, int(counts[a]), total_count, len(ref), motif]
			seqinfo.append(alt)
			add_modified_seq(alignment, sequence, ref, relative_pos)
	raw_seqs = []
	for s in range(1,len(alignment)):
		modified_seq = alignment[s]
		[pos, ref, sequence, count, total_count, max_reach, motif] = seqinfo[s]
		raw_seq = modified_seq.replace('-','')
		duplicate = " "
		if raw_seq in raw_seqs:
			duplicate = "d"
		else:
			raw_seqs.append(raw_seq)
		raw_seq_padded = raw_seq + '-'*(len(modified_seq)-len(raw_seq))
		relative_pos = pos-smallest_ref_pos
		print " %12d  %3d  %s   %s   %s" % (pos, relative_pos+1, modified_seq, raw_seq_padded, duplicate),

		print "    %7d  %7d    %s    %s" % (count, total_count, snip(ref), snip(sequence)),
		if motif:
			print "  motif : %s  %d times" % (motif, number_of_appearances(motif, raw_seq[relative_pos:])),
			print "  in ref %d times" %  number_of_appearances(motif, original_region[relative_pos:]),
		print

	find_chainable([original_region]+raw_seqs, motifs)


	return

#########################################
def get_region_from_das(assembly, chrom, start, end):
	das_request  = "http://genome.ucsc.edu/cgi-bin/das/%s/" % assembly
	das_request += "dna?segment=chr%s:%s,%s" %(chrom, start, end)

	ret = requests.get(das_request).text.lower().replace("\n","")
	pattern = re.compile("<dna.*?>([\w\s]*)</dna>")
	matches = re.findall(pattern, ret)
	dna = ""
	for m in matches: dna += m.replace(" ","")
	return dna

#########################################
def main():

	db, cursor = connect()
	chrom = "1"
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
		[pos, ref, alt, var_counts, total_count, hotspot_id] = variant
		# from these further remove cases where the variant field is a list of SNPs
		list_of_alts = alt.split(",")
		ref_len = len(ref)
		if ref_len==1 and ("," in alt) and (len(alt)+1)==2*len(list_of_alts): continue
		# remove cases with count 0
		new_alts = []
		new_counts = []
		list_of_counts = var_counts.split(",")
		for i in range(len(list_of_alts)):
			if int(list_of_counts[i])==0: continue
			new_alts.append(list_of_alts[i])
			new_counts.append(list_of_counts[i])
		max_reach = pos + ref_len - 1
		candidates.append([pos, ref, ",".join(new_alts),  ",".join(new_counts), total_count, max_reach]) # I want the positions to remain sorted

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

	print "Done clustering."
	reasonable_clusters = [c for c in clusters if len(c)>=2]
	ct = 0
	for cluster in reasonable_clusters[15:18]:
		ct += 1
		reconstruct_alignment(chrom, cluster)

	print ct
	print
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()

# 8,608,914