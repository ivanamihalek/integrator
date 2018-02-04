#!/usr/bin/python
#
#
#

from integrator_utils.python.mysql     import *
from integrator_utils.python.restfuls  import *
from integrator_utils.python.pdb_tools import *
import integrator_utils.python.pdb_constants as pdbc


#########################################
def main():

	if not len(sys.argv)>1:
		print "Usage: %s <seq_file (fasta)> " % sys.argv[0]
		exit(0)


#########################################
if __name__ == '__main__':
	main()
