#!/usr/bin/python

import re
import sys

####################################
def main():

	if len(sys.argv) < 2:
		print  "usage: %s  <infile>" % sys.argv[0]
		exit(1)
	infname  = sys.argv[1]
	infile   = open(infname,'r')
	instring = infile.read().replace('\n', '')
	raw_pattern  = re.compile(r'<tr.*?>(.+?)</tr>')
	data_pattern = re.compile(r'<td.*?>(.+?)</td>')

	for table_row  in re.findall(raw_pattern, instring):
		print "\t".join([re.sub(r'<.*?>', '', data) for data in re.findall(data_pattern, table_row)])
	return

####################################
if __name__ == '__main__':
	main()

