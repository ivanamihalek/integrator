#!/usr/bin/python


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
	number_of_motifs_2 = len(p2)
	# exceptions: not sure how to handle this
	shorter = min(number_of_motifs,number_of_motifs_2)
	if shorter==1:
		return p1[0][0] == p2[0][0]
	if number_of_motifs != len(p2): return False
	matching_motifs = filter(lambda i: p1[i][0] == p2[i][0], range(number_of_motifs))
	return len(matching_motifs) == number_of_motifs

def to_string(pattern):
	string = ""
	for [motif,number_of_repeats] in pattern:
		string += motif*number_of_repeats
	return string

def prettyprint(pattern):

	return ",".join([p[0]+ ":"+ p[1] for p in pattern])


#######################################
###################################

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


