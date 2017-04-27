#!/usr/bin/python

# getting the  genomic region from UCSC DAS server
# https://github.com/ivanamihalek/progesterone/blob/master/PR_binding_sites/scripts/ruby_modules/httputils.rb
# python http request http://docs.python-requests.org/en/master/

from integrator_utils.mysql import *
from integrator_utils.periodic_mods import *

import requests # http  requests
import re



assembly = "hg19"

#########################################
def to_almt_pos(sequence, nongapped_pos):
	almt_pos  = 0
	nongap_ct = 0
	for character in sequence:
		if nongap_ct >= nongapped_pos: break
		if character!="-": nongap_ct +=1
		almt_pos += 1
	return almt_pos

#########################################
def insert (alignment, s, reference_position, ins):
	# insert ins right after reference_pos
	# for example, insert GG after position 5, that is, at positions 6 and 7
	for i in range(len(ins)): # position in the insert
		seq_pos   = reference_position + 1 + i
		almt_pos  = to_almt_pos (alignment[s], seq_pos)
		if almt_pos<len(alignment[s]) and  alignment[s][almt_pos] == "-":
			new_seq = alignment[s][:almt_pos] + ins[i] + alignment[s][almt_pos+1:]
			alignment[s] = new_seq
		else:
			for s2 in range(len(alignment)):
				new_seq = alignment[s2][:almt_pos]
				if s2==s:
					new_seq += ins[i]
				else:
					new_seq += "-"
				new_seq += alignment[s2][almt_pos:]
				alignment[s2] = new_seq
	return

#########################################
def reconstruct_alignment(chrom, cluster):


	periodic = periodicity_found(cluster)
	if periodic: print "periodic"

	smallest_ref_pos = min ([x[0] for x in cluster])
	upper_bound_ref = max ([x[0] + len(x[1]) for x in cluster])
	refseq = "x" * (upper_bound_ref - smallest_ref_pos)
	for v in cluster:
		[pos, ref, alts, var_counts, total_count, max_reach] = v
		relative_pos = pos - smallest_ref_pos
		new_seq = refseq[:relative_pos] + ref + refseq[relative_pos + len(ref):]
		refseq = new_seq
	print " "*19, get_region_from_das(assembly, chrom, smallest_ref_pos, upper_bound_ref).upper()
	print " %12d  %3d  %s" %  (smallest_ref_pos, 1, refseq)
	print "     -------------------------------------------"
	alignment = []
	seqinfo   = []
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
			alt = [pos, ref, sequence, int(counts[a]), total_count, len(ref), motif]
			seqinfo.append(alt)

			refseq = alignment[0]
			if len(sequence) <= len(ref):
				seq_padded   = sequence + "-"*(len(ref)-len(sequence)) # negative padding == no padding
				almt_pos = to_almt_pos(refseq, relative_pos)
				modified_seq = refseq[:almt_pos] + seq_padded + refseq[almt_pos + len(ref):]
				alignment.append(modified_seq)
			else:
				# append the reference sequence, than put an insert
				alignment.append(refseq)
				s = len(alignment)-1 # the index of our sequence
				# insert sequence[len(ref):] right after relative_pos+len(ref)-1
				# for example, for variant A-->AGG at position 5, ref = A,  len(ref) = 1,
				# sequence[len(ref):] = GG and we are inserting it right after position 5
				insert (alignment,  s, relative_pos+len(ref)-1, sequence[len(ref):])

	for s in range(len(alignment)):
		modified_seq = alignment[s]
		[pos, ref, sequence, count, total_count, max_reach, motif] = seqinfo[s]
		print " %12d  %3d  %s" % (pos, pos-smallest_ref_pos+1, modified_seq),
		print "       %7d  %7d    %s    %s" % (count, total_count, ref[:10], sequence[:10]),
		if motif: print "  motif : %s" % motif,
		print

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
		max_reach = pos +max([ref_len]+ [len(x) for x in list_of_alts]) - 1
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

	print "Done clustering."
	ct = 0
	for cluster in clusters:
		if len(cluster)<2: continue
		print "\n********************"
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