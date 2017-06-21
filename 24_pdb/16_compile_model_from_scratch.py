#!/usr/bin/python
#
#
from integrator_utils.mysql import *
import os, subprocess, shutil
import urllib2
from bs4   import BeautifulSoup

swissmodel_dir = "/databases/swissmodel"
scratch_dir    = "/home/ivana/scratch"
struct         = "/home/ivana/projects/enzyme_modeling/code/struct/struct"
pdb_affine     = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_affine_tfm.pl"
extract_chain  = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_extract_chain.pl"
geom_epitope   = "/home/ivana/pypeworks/integrator/24_pdb/integrator_utils/geom_epitope.pl"
blastp         = "/usr/local/bin/blastp"
pdb_blast_db   = "/databases/pdb/blast/pdb_seqres.fasta" # expected to be formatted using makeblastdb
compiled_model_repository = "/home/ivana/monogenic/public/pdb"

conda          = "/home/ivana/miniconda2/bin"
rdkitenv       = "rdkit-env"
rdkit_runner   = "/home/ivana/pypeworks/integrator/24_pdb/17_rdkit_smiles_comparison.py"


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
def substrate_smiles_from_metacyc(cursor, ec_number):
	switch_to_db(cursor,'blimps_development')
	# search by gene may fail because metacyc does not have human version described
	#qry = "select unique_id from metacyc_genes where common_name='%s'" % gene_symbol
	qry  = "select enzymatic_reaction, rxn_left, rxn_right from metacyc_reactions "
	qry += "where ec_number = 'EC-%s'" % ec_number
	ret = search_db(cursor,qry)
	if not ret:
		print "\t failed finding %s in metacyc_reactions "  % ec_number
		return None

	cofactors = []
	alternative_cofactors = []
	alternative_substrates = []
	substrates = []
	for line in ret:
		enzymatic_reaction, rxn_left, rxn_right = line
		for substrate in rxn_left.replace(" ","").split(";") + rxn_right.replace(" ","").split(";"):

			if substrate!='' and '|' not in substrate and substrate not in substrates: substrates.append(substrate)
		for er in enzymatic_reaction.split (";"):
			qry  = "select cofactors, alternative_cofactors, alternative_substrates "
			qry += "from metacyc_enzrxns where unique_id='%s' " % er.replace(" ","")
			ret2 = search_db(cursor,qry)
			for line2 in ret2:
				cofs, alt_cofs, alt_subs = line2
				# list is reserved word, so lets call this array
				for element, array in [[cofs,cofactors], [alt_cofs, alternative_cofactors], [alt_subs, alternative_substrates]]:
					if element!='' and 'EV-EXP' not in element and element not in array: array.append(element)
	subs_smiles = []
	cofactors_smiles = []
	for name, array in [['substrates', substrates], ['cofactors',cofactors],
	                    ['alt_cofs', alternative_cofactors], ['alt_subs', alternative_substrates]]:
		if len(array) == 0: continue
		if 'subs' in name:
			smiles_list = subs_smiles
		else:
			smiles_list = cofactors_smiles
		for compound_id in array:
			qry = "select smiles from metacyc_compounds where unique_id='%s'" % compound_id
			ret2 = search_db(cursor,qry)
			if not ret2: continue
			for line in ret2:
				smiles_list.append([compound_id, line[0]])

	return subs_smiles, cofactors_smiles

##########################################
def download_and_store_ligands_from_pdb(cursor, pdb_id):
	pdb_request = "https://www.rcsb.org/pdb/rest/ligandInfo?structureId=%s" % pdb_id
	response = urllib2.urlopen(pdb_request)
	html = response.read()
	soup = BeautifulSoup(html, 'html.parser')
	fixed_fields = {'pdb_id':pdb_id}
	switch_to_db(cursor, 'monogenic_development')
	for ligand in soup.find_all('ligand'):
		fixed_fields['chemical_id'] = str(ligand['chemicalid'])
		update_fields = {'ligand_type': str(ligand['type'])}
		update_fields['chemical_name'] = str(ligand.chemicalname.string)
		update_fields['smiles'] =  str(ligand.smiles.string)
		store_or_update(cursor, 'pdb_ligands', fixed_fields, update_fields, verbose=True)
	return

##########################################
def  get_pdb_ligand_info(cursor, pdb_id):
	pdb_id = pdb_id.lower()[:4]
	qry = "select  chemical_id, smiles from monogenic_development.pdb_ligands where pdb_id='%s'" % pdb_id
	ret = search_db(cursor,qry)
	if not ret:
		download_and_store_ligands_from_pdb(cursor, pdb_id)
		ret = search_db(cursor,qry)
	if not ret or len(ret)==0:
		return None
	smiles = {}
	for line in ret:
		chemical_id, smile = line
		smiles[chemical_id] = smile
	return smiles

##########################################
def find_pdb_ids_of_similar_seqs(cursor,uniprot_id,scratch):
	os.chdir(scratch)
	qry = "select sequence from monogenic_development.uniprot_seqs where uniprot_id='%s'" % uniprot_id
	ret = search_db(cursor,qry)
	queryfile = uniprot_id+".fasta"
	outf = open(queryfile,"w")
	outf.write(">{}\n{}\n".format(uniprot_id,ret[0][0]))
	outf.close()

	# outfmt 6 is tabular format; the column headers are given using outfmt 7:
	# query acc., subject acc., % identity, alignment length, mismatches,
	#      gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	target_ids = []
	cmd = "{} -db {} -query {} -outfmt 6".format(blastp, pdb_blast_db, queryfile)
	for line in subprocess.check_output(cmd, shell=True).split("\n"):
		field = line.split()
		if len(field)==0 or field[0] != uniprot_id:  continue # this is not the result line
		pct_identity = field[2]
		if float(pct_identity)<40: break
		target_id = field[1]
		if not target_id in target_ids:  target_ids.append(target_id)
	return target_ids

##########################################
# I'll probabily need to optimza this some time soon
def rdkit_compare(smiles1, smiles2):
	cmd = "source {}/activate {}; {} '{}' '{}'; source {}/deactivate {}".format(conda,rdkitenv,rdkit_runner, smiles1, smiles2, conda, rdkitenv)
	# unless specified, the subprocess shell defaults to /bin/sh
	# /bin/sh does not know what is 'source', /bin/bash does
	similarity = 0
	ret = subprocess.check_output(cmd, shell=True, executable="/bin/bash")
	if ret and ret!="": similarity = float(ret.rstrip())
	print smiles1, smiles2, similarity
	return similarity

##########################################
def select_pdbs_with_relevant_ligands(cursor, pdb_id_list, substrate_smiles, cofactor_smiles):
	for pdb_id in pdb_id_list:
		# check whether we already have this pdb in the database - if not download
		smiles = get_pdb_ligand_info(cursor, pdb_id)
		for compound_id, smiles_string in smiles.iteritems():
			# x[0] is the name of the compound, x[1] is the actual smiles
			max_similarity_with_substrate = max([0] + [rdkit_compare(smiles_string,x[1]) for x in substrate_smiles])
			max_similarity_with_cofactor  = max([0] + [rdkit_compare(smiles_string,x[1]) for x in cofactor_smiles])
			if max_similarity_with_substrate>0.5:
				print "possible substrate: %s %s  %5.2f" % (pdb_id, compound_id, max_similarity_with_substrate)
			if max_similarity_with_cofactor>0.5:
				print "possible substrate: %s %s  %5.2f" % (pdb_id, compound_id, max_similarity_with_cofactor)
		#if not soup or not soup.smiles or not soup.dna.string: return None
		#stripped = soup.dna.string.strip().upper()
	exit()



##########################################
def main():
	swissmodel_dir = "/databases/swissmodel"
	db, cursor = connect()
	for dependency in [swissmodel_dir, scratch_dir, struct, pdb_affine,
	                   geom_epitope, extract_chain, compiled_model_repository,
					   blastp, pdb_blast_db, conda, rdkit_runner]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	# first lets focus on the proteins from the newborn screening
	nbs_genes = newborn_screening_genees(cursor)
	for disease in nbs_genes:
		print disease
		for [gene_symbol, ensembl_gene_id, uniprot_id, ec_number] in nbs_genes[disease]:
			print "\t", gene_symbol, ensembl_gene_id, uniprot_id, ec_number
			# check swiss model structure exists, is nonempty, and is indeed pdb
			swissmodel = check_pdb_exists(swissmodel_dir, gene_symbol)
			if not swissmodel: continue
			print "\t", swissmodel
			# is the residue numbering correct?
			identical_pct = check_res_numbers(cursor, "/".join([swissmodel_dir, gene_symbol[0], gene_symbol]), swissmodel)
			mismatch = False
			for chain,pct in identical_pct.iteritems():
				if pct<90:
					print "\t structure seq mismatch? chain", chain, "pct identical positions:",  pct
					mismatch = True
					break
			if mismatch: continue
			# if we got ot here, we might need some scratch space
			scratch = "/".join([scratch_dir,uniprot_id])
			shutil.rmtree(scratch,ignore_errors=True)
			os.makedirs(scratch)

			# what is the biological assembly for this protein? - lets take for now that swissprot is right
			# is this an enzyme?
			if ec_number: # if yes, find substrates and cofactors (ions?)
				# get substrates as smiles strings
				subs_smiles, cofactors_smiles = substrate_smiles_from_metacyc(cursor,ec_number)
				# find all other pdb files with similar sequence and ligands
				pdb_id_list = find_pdb_ids_of_similar_seqs(cursor,uniprot_id,scratch)
				# fid smiles for the ligand - is it similar to the substrate? if not, drop the whole pdb
				# find ligands that are similar to substrates/cofactors, that are not already present in the model
				select_pdbs_with_relevant_ligands(cursor, pdb_id_list, subs_smiles, cofactors_smiles)
			# map ligands onto the model structure - rename them by the nearest chain, and remove if not close
			# map_ligands()
			# store separately distances to substtrates, cofactors, and interfaces
			# shutil.rmtree(scratch,ignore_errors=True)
	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
