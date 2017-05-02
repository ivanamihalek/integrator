#!/usr/bin/python
from collections import deque

#######################################
def decomposition(string, motif):
	# decomposition is left_pad, multiplier (how many times does the motif repeat), right_pad
	# the trivial decomposition, if there is no match inside
	decompositions = [[string, 1, ""]]
	l = len(motif)
	match_positions = deque()
	for i in range(len(string)):
		if string[i:i+l]==motif:
			match_positions.append(i)
	multiplier = 0
	start      = -1
	while match_positions:
		mp = match_positions.popleft()
                if multiplier==0: start = mp
                multiplier += 1
		if len(match_positions)==0 or mp+l != match_positions[0]:
                   decompositions.append([string[:start], multiplier, string[start+multiplier*l:] ])
                   multiplier = 0
                
	return decompositions

#########################################
def main():
    for string in ["lllabcrrr", "lllabcabcrrr" , "lllabc", "abcrrr" ,  "lllrrrrrr", "",
                   "abc", "abcabc",   "lllabcabcabcrrr", "lllabcabcmmmmabcrrr" ]:
        print string
        print decomposition(string, "abc")
        print
    return

#########################################
if __name__ == '__main__':
	main()

