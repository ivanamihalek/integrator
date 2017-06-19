#!/usr/bin/python

# not quite from scratch - use swissmodel so I don't have to model in the loops
# but look for biologically active structure on my own

from integrator_utils.mysql import *
import os, subprocess

swissmodel_dir = "/databases/swissmodel"
scratch_dir    = "/home/ivana/scratch"
struct         = "/home/ivana/projects/enzyme_modeling/code/struct/struct"
pdb_affine     = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_affine_tfm.pl"
extract_chain  = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_extract_chain.pl"
geom_epitope   = "/home/ivana/pypeworks/integrator/24_pdb/integrator_utils/geom_epitope.pl"
compiled_model_repository = "/home/ivana/monogenic/public/pdb"
cwd    = os.getcwd()
chain_pos = 21  #length 1
res_name_pos = 17
res_name_length = 3
res_number_pos = 22
res_number_length = 4
aa_translation = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
# http://www-structmed.cimr.cam.ac.uk/Course/Crystals/optimization.html
# popular crystallographic additives
# * Glycerol, which may stop nucleation and may give you fewer, larger crystals, and
# has the advantage of doubling as a cryo-protectant. "GOL","PGO","PGR"
# * Ethanol or dioxane, which have the effect of poisoning the crystals and stopping too much nucleation
# *Divalent cations like magnesium "EDO","EOH", "DIO
# * A detergent such as beta-octyl-glucoside "SOG","HTG"
# NAG?  N-ACETYL-D-GLUCOSAMINE is that n additive
# http://www.chem.gla.ac.uk/research/groups/protein/mirror/stura/cryst/add.html
crystallographic_additives = ['HOH', "SO4", "GOL","PGO","PGR","EDO","EOH","DIO","SOG","HTG","CL","PEG","NAG", "SAM", "HEM","PE4","ACT","NA"]
physiological_ions =["FE","FE2","MN","ZN","ZN2","MG","CU","CO","CD","MO","VA","NI","W", "SE","CA"]

##########################################
def newborn_screening_genees(cursor):
	nbs_genes = {}
	switch_to_db(cursor,'monogenic_development')
	qry = "select name_short, omim_ids from diseases"
	for line in search_db(cursor, qry):
		name_short, omim_ids = line
		nbs_genes[name_short] = []
		for omim_id in omim_ids.split(';'):
			qry  = "select approved_symbol, ensembl_gene_id from blimps_development.omim_genemaps "
			qry += "where mim_number='%s'" % omim_id
			ret = search_db(cursor,qry, verbose=True)
			approved_symbol, ensembl_gene_id = ret[0]
			qry  = "select uniprot_id, ec_number from blimps_development.uniprot_basic_infos where ensembl_gene_id='%s'" % ensembl_gene_id
			ret = search_db(cursor,qry, verbose=True)
			uniprot_id, ec_number = ret[0]
			nbs_genes[name_short].append([approved_symbol, ensembl_gene_id, uniprot_id, ec_number])
	return nbs_genes

##########################################
def check_pdb_exists(swissmodel_dir, gene_symbol):
	directory = "/".join([swissmodel_dir, gene_symbol[0], gene_symbol])
	if not os.path.exists(directory):
		print "\t", directory," not found"
		return False
	if not os.path.isdir(directory):
		print "\t",directory," not directory"
		return False
	swiss = None
	for model in next(os.walk(directory))[2]:
		if not 'swissmodel' in model: continue
		swiss =  model
	if not swiss:
		print "\t swissmodel not found"
		return False
	cmd = "head -n1 {}".format("/".join([directory,swiss]))
	first_field = subprocess.check_output(cmd, shell=True).split()[0]
	if first_field in ['ATOM','TITLE', 'HEADER','HETATM']:
		return swiss
	print "\t", swiss,"not a PDB file?"
	return False

##########################################
def pdb_to_sequence(path, filename):
	sequence = {}
	infile = open (path+"/"+filename, "r")
	for line in infile:
		if line[:4] != 'ATOM': continue
		chain = line[chain_pos]
		resname = line[res_name_pos:res_name_pos+res_name_length]
		resnumber = int(line[res_number_pos:res_number_pos+res_number_length])
		if not aa_translation.has_key(resname): continue
		if not sequence.has_key(chain): sequence[chain] = {}
		sequence[chain][resnumber] = aa_translation[resname]

	return sequence

##########################################
def check_res_numbers(cursor, path, model):

	field = model.split("_")
	uniprot = field[0]
	qry = "select sequence from monogenic_development.uniprot_seqs where uniprot_id='%s'" % uniprot
	ret = search_db(cursor,qry)
	if not ret:
		print model
		print 'sequence not found for', uniprot
		exit()
	uniprot_sequence = ['.']+list(ret[0][0]) # start index from 1
	unilength = len(uniprot_sequence)
	pdb_sequence    = pdb_to_sequence(path, model)
	identical_count = {}
	identical_pct = {}
	for chain in pdb_sequence:
		identical_count[chain] = 0
		#print chain
		#print pdb_sequence[chain]
		for resnumber, restype in pdb_sequence[chain].iteritems():
			if resnumber<unilength and restype==uniprot_sequence[resnumber]:
				identical_count[chain] += 1
		identical_pct[chain] = float(identical_count[chain])/min(unilength, len(pdb_sequence[chain]))
		identical_pct[chain] = int(100*identical_pct[chain])
		#print identical_pct[chain]
		#print uniprot_sequence
		#exit()
	return identical_pct

##########################################
def main():
	swissmodel_dir = "/databases/swissmodel"
	db, cursor = connect()
	for dependency in [swissmodel_dir, scratch_dir, struct, pdb_affine,
	                   geom_epitope, extract_chain, compiled_model_repository]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	# first lets focus on the proteins from the newborn screening
	nbs_genes = newborn_screening_genees(cursor)
	for disease in nbs_genes:
		print disease
		for [gene_symbol, ensembl_gene_id, uniprot_id, ec_number] in nbs_genes[disease]:
			if gene_symbol!='MUT': continue
			print "\t", gene_symbol, ensembl_gene_id, uniprot_id, ec_number
			# check swiss model structure exists, is nonempty, and is indeed pdb
			swissmodel = check_pdb_exists(swissmodel_dir, gene_symbol)
			if not swissmodel: continue
			print "\t", swissmodel
			# is the residue numbering correct?
			identical_pct = check_res_numbers(cursor, "/".join([swissmodel_dir, gene_symbol[0], gene_symbol]), swissmodel)
			mismatch = False
			for chain,pct in identical_pct.iterites():
				if pct<90:
					print "\t structure seq mismatch? pct ideticap positions:",  pct
					mismatch = True
					break
			if mismatch: continue
			# what is the biological assembly for this protein? - lets take for now that swissprot is right
			# is this an enzyme?
			if ec_number: # if yes, find substrates and cofactors (ions?)
				
			# get substrates as smiles strings
			# find all other pdb files with similar sequence and ligand
			# fid smiles for the ligand - is it similar to the substrate? if not, drop the whole pdb
			# store separately distances to lignds, ions, and interfaces
	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
