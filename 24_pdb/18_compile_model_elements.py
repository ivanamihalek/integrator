#!/usr/bin/python
#
#
from integrator_utils.mysql import *
from integrator_utils.restfuls import *
import os, subprocess, shutil
from collections import OrderedDict

swissmodel_dir = "/databases/swissmodel"
scratch_dir    = "/home/ivana/scratch"
struct         = "/home/ivana/code/struct/struct"
pdb_affine     = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_affine_tfm.pl"
extract_chain  = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_extract_chain.pl"
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
physiological_ions = ["FE","FE2","MN","ZN","ZN2","MG","CU","CO","CD","MO","VA","NI","W", "SE","CA","K"]

##########################################
def newborn_screening_genes(cursor):
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
			qry  = "select uniprot_id, ec_number, cofactors from blimps_development.uniprot_basic_infos where ensembl_gene_id='%s'" % ensembl_gene_id
			ret = search_db(cursor,qry, verbose=True)
			uniprot_id, ec_number,uniprot_cofactors = ret[0]
			nbs_genes[name_short].append([approved_symbol, ensembl_gene_id, uniprot_id, ec_number, uniprot_cofactors])
	return nbs_genes

##########################################
def all_iem_related_genes(cursor):
	iem_genes = {}
	qry  = "select mim_number, approved_symbol, ensembl_gene_id, phenotypes "
	qry += "from blimps_development.omim_genemaps "
	qry += "where inborn_error_of_metabolism=1"
	ret = search_db(cursor,qry, verbose=False)
	if not ret:
		print "no iems ?!"
		exit()
	for row in ret:
		[mim_number, approved_symbol, ensembl_gene_id, phenotypes] = row
		qry  = "select uniprot_id, ec_number, cofactors from blimps_development.uniprot_basic_infos where ensembl_gene_id='%s'" % ensembl_gene_id
		ret = search_db(cursor,qry, verbose=False)
		if not ret:
			print "no uniprot info found for", ensembl_gene_id
			continue
			#exit()
		uniprot_id, ec_number, cofactors = ret[0]
		#print mim_number, approved_symbol, ensembl_gene_id, uniprot_id, ec_number, phenotypes
		if not iem_genes.has_key(phenotypes): iem_genes[phenotypes] = []
		iem_genes[phenotypes].append([approved_symbol, ensembl_gene_id, uniprot_id, ec_number, cofactors])
	return iem_genes

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
	ivana = None
	for model in next(os.walk(directory))[2]:
		if 'swissmodel' in model:  swiss = model
		if 'ivana' in model:  ivana = model
	if not ivana and not swiss:
		print "\t swissmodel not found"
		return False
	if ivana: swiss = ivana # use model made by hand if available
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
		for resnumber, restype in pdb_sequence[chain].iteritems():
			if resnumber<unilength and restype==uniprot_sequence[resnumber]:
				identical_count[chain] += 1
		identical_pct[chain] = float(identical_count[chain])/min(unilength, len(pdb_sequence[chain]))
		identical_pct[chain] = int(100*identical_pct[chain])
	return identical_pct

##########################################
def get_regulators(cursor, reg_by):
	regulators = []
	for reg_id in reg_by.split(";"):
		qry  = "select  regulator from metacyc_regulations where unique_id='%s'" % reg_id
		ret = search_db(cursor,qry)
		if not ret: continue
		regulators.append(ret[0][0])
	return  ";".join(regulators)

##########################################
def substrate_smiles_from_metacyc(cursor, ec_numbers):
	switch_to_db(cursor,'blimps_development')
	rets = []
	# search by gene may fail because metacyc does not have human version described
	#qry = "select unique_id from metacyc_genes where common_name='%s'" % gene_symbol
	for ec_number in ec_numbers.split(";"):
		qry  = "select enzymatic_reaction, rxn_left, rxn_right from metacyc_reactions "
		qry += "where ec_number = 'EC-%s'" % ec_number
		ret = search_db(cursor,qry)
		if ret: rets += ret
	if len(rets)==0:
		print "\t failed finding %s in metacyc_reactions "  % ec_numbers
		return None

	cofactors = []
	alternative_cofactors = []
	alternative_substrates = []
	regulators = []
	substrates = []

	for line in ret:
		enzymatic_reaction, rxn_left, rxn_right = line
		for substrate in rxn_left.replace(" ","").split(";") + rxn_right.replace(" ","").split(";"):
			if substrate in ['WATER', 'NA','CL','PROTON']: continue
			if substrate!='' and '|' not in substrate and substrate not in substrates: substrates.append(substrate)

		for er in enzymatic_reaction.replace(" ","").split(";"):
			if len(er)==0: continue
			qry  = "select cofactors, alternative_cofactors, alternative_substrates, regulated_by "
			qry += "from metacyc_enzrxns where unique_id='%s' " % er.replace(" ","")
			ret2 = search_db(cursor,qry)
			for line2 in ret2:
				cofs, alt_cofs, alt_subs, reg_by = line2
				# regulators need an extra resolution step
				regs = get_regulators(cursor, reg_by)
				# list is reserved word, so lets call this array
				for element, array in [[cofs,cofactors], [alt_cofs, alternative_cofactors],
				                       [alt_subs, alternative_substrates], [regs, regulators]]:

					if element in ['WATER', 'NA','CL','PROTON']: continue
					if element!='' and 'EV-EXP' not in element and element not in array: array.append(element)
	subs_smiles = []
	cofactors_smiles = []
	regs_smiles = []
	smiles_list = {'substrates': subs_smiles,'alt_subs': subs_smiles,
	               'cofactors': cofactors_smiles,'alt_cofs': cofactors_smiles, 'regulators': regs_smiles}
	for name, array in [['substrates', substrates], ['cofactors',cofactors],['alt_cofs', alternative_cofactors],
	                    ['alt_subs', alternative_substrates], ['regulators', regulators]]:
		if len(array) == 0: continue
		for compound_id in array:
			qry  = "select smiles from metacyc_compounds where unique_id='%s'" % compound_id
			ret2 = search_db(cursor,qry)
			if not ret2:
				for comp_id in compound_id.split(";"):
					comp_id = comp_id.replace(" ","").replace("+","").upper()
					smiles_list[name].append([comp_id, None]) # single atoms or ions may not have the smiles string
			else:
				for line in ret2:
					smiles_list[name].append([compound_id, line[0]])

	del smiles_list['alt_subs']
	del smiles_list['alt_cofs']
	return smiles_list

##########################################
def parse_ions(uniprot_cofactors, metacyc_cofactors):
	for cof in uniprot_cofactors.split(";"):
		cof =  cof[:2].upper()
		# TODO: could it be that i have some cofactors here other than ions, which are not present in metacyc?
		if cof not in physiological_ions: continue
		seen_in_metacyc = False
		for [compound, smiles] in metacyc_cofactors:
			if cof==compound[:2].upper():
				seen_in_metacyc = True
				break
		if not seen_in_metacyc: metacyc_cofactors.append([cof,cof])

##########################################
def download_and_store_ligands_from_pdb(cursor, pdb_id):
	return_list = ligands_from_pdb(pdb_id) # restful request
	fixed_fields = {'pdb_id':pdb_id}
	switch_to_db(cursor, 'monogenic_development')
	for ligand_info in return_list:
		fixed_fields['chemical_id'] = ligand_info['chemical_id']
		update_fields = {'ligand_type': ligand_info['ligand_type']}
		update_fields['chemical_name'] = ligand_info['chemical_name']
		update_fields['smiles'] = ligand_info['smiles']
		store_or_update(cursor, 'pdb_ligands', fixed_fields, update_fields, verbose=False)
	return

##########################################
def get_pdb_ligand_info(cursor, pdb_id):
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
	target_pct_idtty = OrderedDict()
	cmd = "{} -db {} -query {} -outfmt 6".format(blastp, pdb_blast_db, queryfile)
	for line in subprocess.check_output(cmd, shell=True).split("\n"):
		print line
		field = line.split()
		if len(field)==0 or field[0] != uniprot_id:  continue # this is not the result line
		pct_identity = field[2]
		if float(pct_identity)<40: break
		target_id = field[1].replace("_","")
		if target_id[-1]==target_id[-2]: target_id=target_id[:-1]
		if not target_id in target_pct_idtty.keys(): target_pct_idtty[target_id]=pct_identity
	return target_pct_idtty

##########################################
def rdkit_compare(smiles1, smiles2):
	cmd = "source {}/activate {}; {} '{}' '{}'; source {}/deactivate {}".format(conda,rdkitenv,rdkit_runner, smiles1, smiles2, conda, rdkitenv)
	# unless specified, the subprocess shell defaults to /bin/sh
	# /bin/sh does not know what is 'source', /bin/bash does
	similarity = 0
	ret = subprocess.check_output(cmd, shell=True, executable="/bin/bash")
	if ret and ret!="": similarity = float(ret.rstrip())
	return similarity

##########################################
def reorganize_ligands(smiles_dict):
	# some ligands might have multiple roles
	# in particular they might be a substrate and a regulator (witness PHE in PAH)
	# on the other hand, there might be some notational confusion regarding what is
	# a substrate and wht is a cofactor
	ligands_reorganized = {}
	for [ligand, smiles] in smiles_dict['substrates']:
		ligands_reorganized[ligand] = [smiles, 'substrate']

	for [ligand, smiles] in smiles_dict['regulators']:
		if ligands_reorganized.has_key(ligand):
			ligands_reorganized[ligand] = [smiles, 'substrate;regulator']
		else:
			ligands_reorganized[ligand] = [smiles, 'regulator']

	for [ligand, smiles] in smiles_dict['cofactors']:
		if ligands_reorganized.has_key(ligand):
			[smiles0, function0] = ligands_reorganized[ligand]
			if function0=='substrate': continue
			ligands_reorganized[ligand] = [smiles, 'cofactor;regulator']
		ligands_reorganized[ligand] = [smiles, 'cofactor']

	return ligands_reorganized

##########################################
def get_ec_from_pdb(cursor, pdb_chain_id):
	switch_to_db(cursor,'monogenic_development')
	qry = "select  uniprot_id, ec_numbers, cofactors from pdb_uniprot_ec_maps where pdb_chain='%s'" % pdb_chain_id
	ret = search_db(cursor,qry)
	if ret: return ret[0] # uniprot, ec numbers, cofactors
	# fetch uniprot id from PDB's das service
	uniprot_id = uniprot_from_pdb_chain(pdb_chain_id)  # restful request
	# fetch enzyme info from Uniprot
	ec_numbers, cofactors = ec_cofactors_from_uniprot(uniprot_id)
	fixed_fields = {'pdb_chain':pdb_chain_id}
	update_fields = {'uniprot_id':uniprot_id, 'ec_numbers': ec_numbers, 'cofactors':cofactors}
	store_or_update(cursor, 'pdb_uniprot_ec_maps', fixed_fields, update_fields)
	return  uniprot_id, ec_numbers, cofactors

##########################################
def get_native_ligands(cursor, pdb_chain_id):
	# check for pdb ec map in the database
	# if not present, fetch from PDB and uniprot
	uniprot_id, ec_numbers, uniprot_cofactors = get_ec_from_pdb(cursor, pdb_chain_id)
	smiles_dict = substrate_smiles_from_metacyc(cursor,ec_numbers)
	if not smiles_dict: smiles_dict = {'cofactors':""}
	parse_ions(uniprot_cofactors,smiles_dict['cofactors'])
	if not smiles_dict: return
	native_ligands_smiles = reorganize_ligands(smiles_dict)
	return native_ligands_smiles

##########################################
def select_pdbs_with_relevant_ligands(cursor, pdb_id_list):

	usable_pdbs = {}
	for pdb_chain_id in pdb_id_list:
		print "\t\t", pdb_chain_id
		# the ligands that this enzyme  processes
		native_ligands_smiles = get_native_ligands(cursor, pdb_chain_id)
		if not native_ligands_smiles: continue
		# the ligands that are actually present in the PDB
		pdb_ligand_smiles = get_pdb_ligand_info(cursor, pdb_chain_id)
		if not pdb_ligand_smiles: continue

		ligand_similarities = []
		for pdb_compound_id, smiles_string in pdb_ligand_smiles.iteritems():
			if pdb_compound_id in crystallographic_additives: continue
			# function: substrate, cofactor, regulator
			for native_ligand_id, [native_smiles, native_ligand_function] in native_ligands_smiles.iteritems():
				# are the two strings trivially equal by any chance?
				similarity=0
				if pdb_compound_id in physiological_ions:
					if not native_ligand_id in physiological_ions: continue
					# pdb people can be lax with the ion charge
					if (pdb_compound_id.translate(None,"+-123"))==(native_ligand_id.translate(None,"+-123")):
						similarity=1.0
				else:
					if native_ligand_id in physiological_ions: continue
					if native_smiles and smiles_string:
						if native_smiles.upper()==smiles_string.upper():
							similarity=1.0
						else:
							similarity = rdkit_compare(native_smiles, smiles_string)
				if similarity>0.5:
					# get rid of the plural s in in the native function
					ligand_similarities.append([pdb_compound_id, native_ligand_id, native_ligand_function, "%.2f"%similarity])
				print "\t\t", pdb_compound_id, native_ligand_id, native_ligand_function, "%.2f"%similarity
			# there can be ligands with the same chemical structure at multiple places in the  sturcture - sort this out when
			# mapping on the model structure
		if len(ligand_similarities)>0:
			usable_pdbs[pdb_chain_id]=ligand_similarities

	return usable_pdbs

##########################################
def print_smiles (name, smiles):
	print "\t%s smiles from metacyc:" %name

	print "\t"+"\n\t".join(map(lambda x: "%s  %s"%(x[0],x[1]),smiles))+"\n"

##########################################
def process_enzyme(cursor,gene_symbol, ensembl_gene_id, ec_number, swissmodel, uniprot_cofactors,uniprot_id,scratch):
	# find all other pdb files with similar sequence and ligands
	print "\tlooking for pdb with sequence similar to  %s ..." % gene_symbol
	pdb_pct_similarity = find_pdb_ids_of_similar_seqs(cursor,uniprot_id,scratch)

	print "\tsearching for pdbs with ligands similar to one of the substrates from metacyc"
	usable_pdbs = select_pdbs_with_relevant_ligands(cursor,pdb_pct_similarity.keys())
	if not usable_pdbs: return
	for usable_pdb_id, sims in usable_pdbs.iteritems():
		print "\t\t", usable_pdb_id, pdb_pct_similarity[usable_pdb_id], sims
		for sim in sims:
			[pdb_ligand_name, metcyc_ligand_name, ligand_function, tanimoto] = sim
			print pdb_ligand_name, metcyc_ligand_name, tanimoto
			fixed_fields  = {'gene_symbol':gene_symbol, 'pdb_chain':usable_pdb_id,
			                 'pdb_ligand':pdb_ligand_name, 'metacyc_ligand':metcyc_ligand_name}
			update_fields = {'ensembl_gene_id':ensembl_gene_id,'uniprot_id':uniprot_id,
			                 'ec_number':ec_number,'main_model':swissmodel,
			                 'pdb_pct_identical_to_uniprot': pdb_pct_similarity[usable_pdb_id],
			                 'ligand_function':ligand_function,
			                 'ligand_tanimoto':tanimoto}
			# just store in the database and continue in the next pipe in the pipeline
			store_or_update(cursor, 'monogenic_development.model_elements', fixed_fields, update_fields, verbose=False)

##########################################
def main():
	swissmodel_dir = "/databases/swissmodel"
	db, cursor = connect()
	for dependency in [swissmodel_dir, scratch_dir, struct, pdb_affine,
	                   extract_chain, compiled_model_repository,
					   blastp, pdb_blast_db, conda, rdkit_runner]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()
	# first lets focus on the proteins from the newborn screening
	genes = newborn_screening_genes(cursor)
	#genes = all_iem_related_genes(cursor)
	for disease in genes.keys():
		print disease
		for [gene_symbol, ensembl_gene_id, uniprot_id, ec_number, uniprot_cofactors] in genes[disease]:
			if gene_symbol=='PAH': continue
			print disease
			print "\t", gene_symbol, ensembl_gene_id, uniprot_id, ec_number, uniprot_cofactors
			#qry = "select * from monogenic_development.model_elements where gene_symbol='%s'" % gene_symbol
			#ret = search_db(cursor,qry)
			#if ret:
			#	print "model elements found"
			#	continue
			# check swiss model structure exists, is nonempty, and is indeed pdb
			print "\tchecking model exists"
			swissmodel = check_pdb_exists(swissmodel_dir, gene_symbol)
			if not swissmodel: continue
			print "\tfound swissmodel:", swissmodel
			# is the residue numbering correct?
			print "\tchecking the numbering in the pdb file ..."
			identical_pct = check_res_numbers(cursor, "/".join([swissmodel_dir, gene_symbol[0], gene_symbol]), swissmodel)
			mismatch = False
			for chain,pct in identical_pct.iteritems():
				print "\t", chain, "pct identitiy", pct
				if pct<90:
					print "\t structure seq mismatch? chain", chain, "pct identical positions:",  pct
					mismatch = True
					break
			if mismatch: continue
			print "\tsequence ok; creating scratch directory\n"
			# if we got ot here, we might need some scratch space
			scratch = "/".join([scratch_dir,uniprot_id])
			shutil.rmtree(scratch,ignore_errors=True)
			os.makedirs(scratch)

			# what is the biological assembly for this protein? - lets take for now that swissprot is right
			# is this an enzyme?
			if ec_number: # if yes, find substrates and cofactors (ions?)
				process_enzyme(cursor,gene_symbol, ensembl_gene_id, ec_number, swissmodel, uniprot_cofactors,uniprot_id,scratch)
			# else if TM protein
			# else
			shutil.rmtree(scratch,ignore_errors=True)
	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
