#!/usr/bin/python


##########################################
def main():
	inf = open("/databases/agilent/v4/S03723314_Regions.bed","r")
	for line in inf:
		fields = line.rstrip().split("\t")
		if len(fields)<3: continue
		[chrom, start, end] = fields[:3]
		if not chrom[:3]=='chr': continue
		print "\t".join([chrom[3:], start,end])
	return True

#########################################
if __name__ == '__main__':
	main()
