#!/usr/bin/python

from integrator_utils.periodic_mods import *


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

