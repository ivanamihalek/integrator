#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess
from math import log

alignment_repository = "/home/ivana/monogenic/public/alignments"

def column_entropy(string):
	bin = {}
	for char in string:
		if not bin.has_key(char): bin[char]=0
		bin[char] +=1
	entropy = 0.0
	for char, count in bin.iteritems():
		freq = float(count)/len(string)
		if count>0:
			entropy -= freq*log(freq)
	return entropy/log(20)

##########################################
def main():

	for dependency in [alignment_repository]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()
	db, cursor = connect()
	qry = "select * from monogenic_development.diseases"
	ret = search_db(cursor, qry)
	for line in ret:
		[id, name_short, name_long, omim_ids, description, prim, sec] = line
		print "="*60
		print name_short
		for omim_id in omim_ids.split(";"):
			qry = "select approved_symbol, ensembl_gene_id from omim_genemaps where mim_number='%s'" % omim_id
			ret2 = search_db(cursor, qry)
			for line2 in ret2:
				approved_symbol = line2[0]
				ensembl_gene_id = line2[1]
				qry = "select uniprot_id from uniprot_basic_infos where gene_name='%s'" % approved_symbol
				ret3 = search_db(cursor, qry)
				if not ret3:
					print "uniprot id not found"
					exit()
				uniprot_id = ret3[0][0]
				infile = open("{}/{}.restr.afa".format(alignment_repository,uniprot_id),"r")
				seq = {}
				for line in  infile:
					line = line.rstrip()
					if len(line)==0: continue
					if line[0]=='>':
						seqname = line[1:]
					else:
						if not seq.has_key(seqname): seq[seqname] =""
						seq[seqname] += line
				infile.close()
				del seq['human']
				columns = ['']*len(seq[uniprot_id])
				for name, sequence in seq.iteritems():
					for i in range(len(sequence)): columns[i]+=sequence[i]
				conservation = []
				for column in columns:
					conservation.append(1-column_entropy(column))
				conservation_score = []
				for i in range(len(columns)):
					column = columns[i]
					start = max(i-2,0)
					end  = min(i+3,len(columns))
					avg = reduce (lambda x,y: x+y, conservation[start:end])/len(conservation[start:end])
					aa_types_present = "".join(set(list(column.replace("-",""))))
					conservation_score.append( "%d:%s:%.3f" % (i+1, aa_types_present, conservation[i] + 0.1*avg) )

				cons_string = ",".join(conservation_score)
				qry  = "update monogenic_development.uniprot_seqs "
				qry += "set conservation_in_vertebrates='%s' " % cons_string
				qry += "where uniprot_id='%s' " % uniprot_id
				search_db(cursor, qry, verbose=False)



########################################
if __name__ == '__main__':
	main()
