#!/usr/bin/python
from collections import deque


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

#########################################
def main():
	if True:
		for string in ["lllabcrrr", "lllabcabcrrr" , "lllabc", "abcrrr" ,  "lllrrrrrr", "",
					   "abc", "abcabc",   "lllabcabcabcrrr", "lllabcabcmmmmabcrrr" ]:
			print string
			print decomposition(string, ["abc"])
			print

		for string in [u'CCAGTCTTT', u'CGTCTTT', u'CCA', u'CCAGTCT']:
			print string
			print decomposition(string, ["T"])
			print

		for string in ['TTCA','TTGGGCA']:
			print string
			print decomposition(string, ["G"])
			print

	for string in ['TGT','TT','TG','TGTT','TTT']:
		print string
		print decomposition(string, ["TG","T"])
		print


	for string in ['CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGC',
				   'CTGCTGTTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGC',
				   'CTGCTGTTGCTGTTGCTGTTGCTGCTGCTGCTGCTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGC']:
		print string
		print decomposition(string, ["TGCTGT", "TGC"])
		print


	return

#########################################
if __name__ == '__main__':
	main()

