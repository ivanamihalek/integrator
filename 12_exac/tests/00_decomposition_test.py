#!/usr/bin/python
from collections import deque

#######################################
def decomposition(string, motif):

        # the trivial decomposition, if there is no match inside
        decomps = [[string, 0, ""]]
        
	# decomposition is left_pad, multiplier (how many times does the motif repeat), right_pad
	l = len(motif)
	match_positions = deque()
	for i in range(len(string)):
		if string[i:i+l]==motif:
			match_positions.append(i)

        if len(match_positions)==0:
               return decomps

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

#########################################
def main():

     if False:
            for string in ["lllabcrrr", "lllabcabcrrr" , "lllabc", "abcrrr" ,  "lllrrrrrr", "",
                           "abc", "abcabc",   "lllabcabcabcrrr", "lllabcabcmmmmabcrrr" ]:
                print string
                print decomposition(string, "abc")
                print

            for string in [u'CCAGTCTTT', u'CGTCTTT', u'CCA', u'CCAGTCT']:
                print string
                print decomposition(string, "T")
                print

            for string in ['TTCA','TTGGGCA']:
                print string
                print decomposition(string, "G")
                print
        
     for string in ['TGT','TT','TG','TGTT','TTT']:
        print string
        print decomposition(string, "T")
        print

        
     return

#########################################
if __name__ == '__main__':
	main()

