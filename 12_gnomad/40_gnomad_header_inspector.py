#!/usr/bin/python3

import shlex


##########################################
def main():
	infile = open("/storage/databases/gnomad/sample_header.txt")
	for line in infile:
		if line[:len("##reference")] == "##reference":
			print("reference assembly:",line.rstrip().split("=")[-1])
			continue
		if line[:6] == '#CHROM': break
		if line[:len('##INFO=<')] != '##INFO=<': continue
		my_splitter = shlex.shlex(line.rstrip()[len('##INFO=<'):-1], posix=True)
		my_splitter.whitespace += ','
		my_splitter.whitespace_split = True
		fields = list(my_splitter)
		field_name =  fields[0].split('=')[-1]
		field_content = fields[-1].split('=')[-1]
		print(field_name)
		# if field_name=='CSQ': # it looks like somehere between versions 2.0.1 and 2.1.1 we gave up that field
		if field_name=='vep': # now we are taking annotation from Ensembl's VEP
			print(field_name,field_content.split(':')[0])
			for subcontent in field_content.split(':')[1].split('|'):
				print("\t", subcontent)
			# however, here, I think we'll do our own annotation
		else:
			print(field_name, "      ",  field_content)

	print()

	return
#########################################
if __name__ == '__main__':
	main()
locals()