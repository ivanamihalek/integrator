#!/usr/bin/python

# getting the  genomic region from UCSC DAS server
# https://github.com/ivanamihalek/progesterone/blob/master/PR_binding_sites/scripts/ruby_modules/httputils.rb
# python http request http://docs.python-requests.org/en/master/

from integrator_utils.mysql import *

import requests # http  requests
import re
from integrator_utils.periodic_mods import *
from time import time

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

#########################################
def characterize_region(cluster):
	smallest_ref_pos = min([x[0] for x in cluster])
	upper_bound_ref  = max([x[0] + len(x[1]) -1 for x in cluster])
	number_of_variants = 0
	for v in cluster:
		[pos, ref, alts, var_counts, total_count, max_reach] = v
		number_of_variants += len( filter(lambda x: len(x)>0 and int(x)>0, var_counts.split(",")) )
	return [smallest_ref_pos, upper_bound_ref, number_of_variants]


#########################################
def get_region_from_das(assembly, chrom, start, end):
	"""
	"""
	das_request  = "http://genome.ucsc.edu/cgi-bin/das/%s/" % assembly
	das_request += "dna?segment=chr%s:%s,%s" %(chrom, start, end)

	ret = requests.get(das_request).text.lower().replace("\n","")
	pattern = re.compile("<dna.*?>([\w\s]*)</dna>")
	matches = re.findall(pattern, ret)
	dna = ""
	for m in matches: dna += m.replace(" ","")
	return dna

 #########################################
def motif_report(chrom, cluster):
	""" find motifs in clusters

	"""
	smallest_ref_pos = min([x[0] for x in cluster])
	upper_bound_ref  = max([x[0] + len(x[1]) - 1 for x in cluster])
	original_region = get_region_from_das(assembly, chrom, smallest_ref_pos, upper_bound_ref).upper()
	refseq = original_region
	if not refseq or len(refseq)==0:
		print " problem obtaining region on the reference sequence:"
		print assembly, chrom, smallest_ref_pos, upper_bound_ref
		return None, None

	seqinfo   = []
	motifs    = []
	alignment = []
	# I don't really need the alignment here, but since I have it implemented ...
	alt = [smallest_ref_pos, refseq[0], refseq[0], 1, 1, 1, None]
	alignment.append(refseq)
	seqinfo.append(alt)
	for v in cluster:
		[pos, ref, alts, var_counts, total_count, max_reach] = v
		if alts=='' or var_counts=='': continue
		relative_pos = pos - smallest_ref_pos
		counts = var_counts.split(",")
		alternatives = alts.split(",")
		a = -1 #  index of the alternative sequence fo this position
		for sequence in alternatives:
			a += 1
			motif = find_motif_in_pair(sequence, ref)
			if motif and not motif in motifs: motifs.append(motif)
			alt = [pos, ref, sequence, int(counts[a]), total_count, len(ref), motif]
			seqinfo.append(alt)
			add_modified_seq(alignment, sequence, ref, relative_pos)
	raw_seqs = [original_region]
	freq = {}
	freq[original_region] = "1:1"
	for s in range(1,len(alignment)):
		modified_seq = alignment[s]
		[pos, ref, sequence, count, total_count, max_reach, motif] = seqinfo[s]
		raw_seq = modified_seq.replace('-','')
		# it seemed to me at one point that eac might have duplicates, but they seem
		# to have caught them and assigned them a frequency of 0
		if count==0: continue
		if raw_seq in raw_seqs:
			pass
		else:
			raw_seqs.append(raw_seq)
		freq[raw_seq] = "%d:%d" % (count,total_count)

	########### REPORT/STORE TO DB ##########################
	# never mind the clustering - I am not sure any more that it helps
	# motif counting though should protect me from counting as diffent
	# variants when a motif indel is assigned to a different place in the repeat expansion
	# as different a variant
	#for cluster in find_comparable_patterns(raw_seqs, motifs):
	report_items= []
	for seq in raw_seqs:
		pattern = decomposition(seq, motifs)
		report_items.append(prettyprint(pattern) + "," + freq[to_string(pattern)])


	return ",".join(motifs), ";".join(report_items)


#########################################
def find_clusters_of_candidates(candidates):
	""" find places with indels that cluster and possibly overlap
		two positions belong to a cluster if they are within 3nt from each other,
		or within each other's reach
		the assumptions is that candidates come ordered by increasing position
		return clusters with more than one member
	"""
	print "no candidates:", len(candidates)
	clusters = []
	for candidate in candidates:
		cluster_found = False
		for cluster in clusters:
			if len([x for x in cluster if x[0]<=candidate[0] <= max(x[0]+3,x[-1])]):
				cluster.append(candidate)
				cluster_found = True
				break
		if not cluster_found: # start new cluster - cluster is a list of candidates
			clusters.append([candidate])

	# lets look at isolated cases too - gnomad might already
	# have some periodic expansions there
	reasonable_clusters = [c for c in clusters if len(c)>=1]
	return reasonable_clusters

########################################
def find_complex_variants(cursor, table):
	""" find locations with frequent indels (not SNVS, but larger)
		I am skipping cases where reference length > 30 (presumably large deletion or indel)
		return a list of candidates, each candindate = [pos, ref, alts, counts, total_count, max_reach]
		also return the number of long indels we have dropped
	"""
	candidates = []
	qry  = "select position, reference, variant, variant_count, total_count from %s " % table
	qry += "where char_length(reference)!=1 or char_length(variant)!=1"
	print qry
	long_vars_count = 0
	for variant in search_db(cursor, qry):
		[pos, ref, alt, var_count, total_count] = variant
		# I don't want large rearrangements here
		if len(ref) > 30:
			long_vars_count += 1
			continue
		# from these further remove cases where the variant field is a list of SNPs
		if alt == '' or var_count == '': continue
		ref_len = len(ref)
		# remove cases with count 0
		max_reach = pos + ref_len - 1
		# I want the positions to remain sorted
		candidates.append([pos, ref, alt, var_count, total_count, max_reach])
	return candidates, long_vars_count


########################################
# the clustering and motif searching bomb for  the region
# 14:106330014-106330121 resulting in the report string > 65,000 chracters
#########################################
def main():
	"""   main
			loop over chroms
			find_complex_variants
			find_clusters_of_candidates
			find periodic motifs in clusters, if any
			create 'report' for each periodic region
			store
	"""
	db, cursor = connect()
	#chroms = ['1','22']
	#chroms = ['2','21']
	#chroms = ['3','20']
	#chroms = ['4','19']
	#chroms = ['5','18']
	#chroms = ['6','17']
	#chroms = ['7','16']
	#chroms = ['8','15']
	#chroms = ['9','14']
	#chroms = ['10','13']
	chroms = ['11','12']
	#chroms = [str(i) for i in range(10,23)]
	#chroms = ['X','Y']
	chroms.reverse()
	for chrom in chroms:
		t0 = time()
		table = "gnomad_freqs_chr_" + chrom
		print
		print "*"*20
		print table
		print "number of variants:", search_db(cursor, "select count(1) from %s" % table)[0][0]
		qry  = "select count(1) from %s " % table
		qry += "where char_length(reference)=1 and char_length(variant)=1"
		print "simple SNPs",  search_db(cursor, qry)[0][0]

		candidates, long_vars_ct = find_complex_variants(cursor, table)
		print
		print "Complex variants with reference<30:", len(candidates),
		print "  long variants: ", long_vars_ct

		clusters =  find_clusters_of_candidates(candidates)
		print
		print "Done clustering. Max pos:", max([cluster[0][0] for cluster in clusters])
		print "Number of hotspot regions:", len(clusters)


		number_of_vars_in_clusters = 0
		number_of_clusters_with_periodic_motifs = 0
		for cluster in clusters:
			# no varaints: cluster is just the number of positions here, not the number of
			# vars repoted for each
			[start,end, number_of_variants] = characterize_region(cluster)
			if number_of_variants<2: continue
			number_of_vars_in_clusters += number_of_variants
			fixed_fields  = {'chrom':chrom, 'start':start, 'end':end}
			store_without_checking(cursor, 'gnomad_hotspots', fixed_fields)
		print
		print "Number of variants with clusters:", number_of_vars_in_clusters
		print "Number of clusters with periodic motifs:", number_of_clusters_with_periodic_motifs
		print
		print "time taken %.2f min" % ((time() - t0) / 60.0)
		print
	cursor.close()
	db.close()

	return
#########################################
if __name__ == '__main__':
	main()
