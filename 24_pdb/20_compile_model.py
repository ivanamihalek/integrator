#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess
import json

swissmodel_dir    = "/databases/swissmodel"
pdb_path          = "/databases/pdb/structures"
scratch_dir       = "/home/ivana/scratch"
struct            = "/home/ivana/code/struct/struct"
pdb_affine        = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_affine_tfm.pl"
pdbdown           = "/home/ivana/pypeworks/integrator/integrator_utils/pdbdownload.pl"
pdb_chain_rename  = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_chain_rename.pl"
extract_chain     = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_extract_chain.pl"
geom_overlap      = "/home/ivana/pypeworks/integrator/integrator_utils/geom_overlap.pl"
geom_epitope      = "/home/ivana/pypeworks/integrator/integrator_utils/geom_epitope.pl"
compiled_model_repository = "/home/ivana/monogenic/public/pdb"
cwd    = os.getcwd()

chain_pos      = 21  #length 1
res_name_pos   = 17
res_name_length = 3
res_number_pos = 22
res_number_length = 4
aa_translation = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                  'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
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
physiological_ions = ["FE","FE2","MN","ZN","ZN2","MG","CU","CO","CD","MO","VA","NI","W", "SE","CA"]

transform = {}

#########################################
def extract_trivial(infile_path, scratch, chain, ligand_resn):
	infile = open (infile_path, "r")
	outstring = ""
	for line in infile:
		if line[:6] != 'HETATM': continue
		if line[chain_pos] != chain:
			continue
		resname = line[res_name_pos:res_name_pos+res_name_length].replace(' ','')
		if resname != ligand_resn: continue
		outstring += line
	infile.close()

	if len(outstring)==0:
		return None
	filename = infile_path.split("/").pop()
	outfilename = "/".join([scratch, filename.replace('.pdb','')+".{}.{}.pdb".format(chain,ligand_resn)])
	outfile = open(outfilename,"w")
	outfile.write(outstring)
	outfile.close()
	return outfilename

#########################################
def extract_peptide (infile_path, chain, chainfile):
	infile = open (infile_path, "r")
	outf   = open (chainfile,"w")
	for line in infile:
		if line[:4]=='ATOM' and line[chain_pos]==chain:
			outf.write(line)
	outf.close()
	infile.close()

#########################################
def extract_ligand_without_chain_name(path, scratch, filename, chain, ligand_resn):
	print "extracting nonames"
	exit()
	# extract all chains
	infile_path = path+"/"+filename
	# extract all ligands
	ligands = extract_ligands(infile_path, chain, ligand_resn)
	print ligands
	exit()
	# extract all chains
	infile  = open (infile_path, "r")
	chains = set()
	for line in infile:
		if line[:4]=='ATOM': chains.add(line[chain_pos])
	infile.close()
	chains = sorted(list(chains))
	for chain in chains:
		chainfile =  "/".join([scratch,"tmp"+chain+".pdb"])
		extract_peptide(infile_path, chain, chainfile)
	# extract ligands
	# assign  ligand to the chain where it leaves the largest footprint
	exit()

#########################################
def extract_ligands(infile_path, scratch,  ligand_resn):
	infile = open (infile_path, "r")
	outfile = {}
	chains = []
	for line in infile:
		if line[:6]!='HETATM': continue
		resname = line[res_name_pos:res_name_pos+res_name_length].replace(' ','')
		if resname!=ligand_resn: continue
		res_number = line[res_number_pos:res_number_pos+res_number_length].replace(' ','')
		if not outfile.has_key(res_number):
			filename = infile_path.split("/").pop()
			chain = line[chain_pos]
			if chain != ' ':
				outf_name= "/".join([scratch, filename.replace('.pdb','')+".{}.{}.pdb".format(res_number,chain)])
				chains.append(chain)
			else:
				outf_name= "/".join([scratch, filename.replace('.pdb','')+".{}.pdb".format(res_number)])
			print outf_name
			outfile[res_number] = open(outf_name,"w")
		outfile[res_number].write(line)

	for outf in outfile.values():outf.close()
	infile.close()

	return [outf.name for outf in outfile.values()],chains

#########################################
def extract_ligand(path, scratch, filename, chain, ligand_resn):
	infile_path = path+"/"+filename
	outfilename = extract_trivial(infile_path, scratch, chain, ligand_resn)
	if outfilename: return outfilename, None
	# otherwise we failed - for example ligand is not associated with the chain we were told to use
	ligands, chains = extract_ligands(infile_path, scratch,  ligand_resn)
	# say we return the first one
	if len(chains): return ligands[0], chains[0]
	exit()
	# still other possibility: ligand has no chainname, ever
	extract_ligand_without_chain_name(path, scratch, filename, chain, ligand_resn)
	return outfilename

#########################################
def map_ligand_to_model(main_model_info,pdb_path,ligand_containing_structure, ligand_chain, ligand_file_path, scratch):

	ligand_file_tfmd_paths = []

	struct_scratch = scratch + "/" + "struct_scratch"
	shutil.rmtree(struct_scratch,ignore_errors=True)
	os.makedirs(struct_scratch)
	os.chdir(struct_scratch)

	# find the  first chain in anchor
	[path, main_model, main_model_chains] = main_model_info
	for main_model_chain in main_model_chains:
		# do we have the tfm already by any chance
		tfm_key = "{}_{}_{}_{}".format(main_model, main_model_chain,ligand_containing_structure, ligand_chain)
		if transform.has_key(tfm_key):
			outf = open("tmp.tfm","w")
			outf.write(transform[tfm_key])
			outf.close()
		else:
			from_structure_path = pdb_path+"/"+ligand_containing_structure
			to_structure_path   = path+"/"+main_model
			cmd = "%s -from %s -c1 %s -to %s -c2 %s" % (struct, from_structure_path, ligand_chain, to_structure_path, main_model_chain)
			subprocess.call(cmd, shell=True)
			almt_file = "{}{}_to_{}{}.0.aln"\
				.format(ligand_containing_structure.replace('.pdb',''), ligand_chain, main_model.replace('.pdb',''),main_model_chain)
			if not os.path.exists(almt_file) or os.stat(almt_file).st_size == 0: return None
			cmd = "grep tfm {} -A 3 | tail -n 3 | sed 's/{}//g' > tmp.tfm".format(almt_file,"%")
			subprocess.call(cmd, shell=True)
			if not os.path.exists('tmp.tfm') or os.stat('tmp.tfm').st_size == 0: return None
			inf = open("tmp.tfm","r")
			transform[tfm_key] = inf.read()
			inf.close()

		ligand_file_tmp_path = ligand_file_path.replace('.pdb','') + ".tmp.pdb"
		cmd = "{} {} tmp.tfm > {}".format(pdb_affine, ligand_file_path, ligand_file_tmp_path)
		subprocess.call(cmd, shell=True)
		# rename chain to match  chain we mapped is to
		ligand_file_tfmd_path = ligand_file_path.replace('.pdb','') + ".tfmd.{}.pdb".format(main_model_chain)
		cmd = "{} {} {} {} > {}".format(pdb_chain_rename, ligand_file_tmp_path, ligand_chain, main_model_chain, ligand_file_tfmd_path)
		subprocess.call(cmd, shell=True)

		ligand_file_tfmd_paths.append(ligand_file_tfmd_path)

	os.chdir(scratch)
	shutil.rmtree(struct_scratch,ignore_errors=True)

	return ligand_file_tfmd_paths

#########################################
def find_last_res_number(pdb_file_path):
	resnumber_atom   = -1
	resnumber_hetatm = -1
	cmd = "awk '$1==\"ATOM\"' {} | tail -n1".format(pdb_file_path)
	ret = subprocess.check_output(cmd, shell=True)
	if ret and len(ret)>res_number_pos+res_number_length: resnumber_atom = int(ret[res_number_pos:res_number_pos+res_number_length])
	cmd = "awk '$1==\"HETATM\"' {} | tail -n1".format(pdb_file_path)
	ret = subprocess.check_output(cmd, shell=True)
	if ret and len(ret)>res_number_pos+res_number_length: resnumber_hetatm = int(ret[res_number_pos:res_number_pos+res_number_length])

	max_resnumber = max(resnumber_atom,resnumber_hetatm)
	return max_resnumber

#########################################
def merge_ligands(path, anchor, compiled_ligand_file_path, ligand_file_tfmd_paths, scratch):

	compiled_ligands_exists = os.path.exists(compiled_ligand_file_path) and os.path.getsize(compiled_ligand_file_path)>0
	os.chdir(scratch)
	# do I have anything closer than1 1 A?
	# note that I am taking the footprint on the ligand, to remove clashing molecules
	if compiled_ligands_exists:
		clashing_ligands = []
		for ligand_file_tfmd_path in ligand_file_tfmd_paths:
			cmd = "{} {} {}".format(geom_overlap, ligand_file_tfmd_path, compiled_ligand_file_path)
			ret = subprocess.check_output(cmd, shell=True).rstrip()
			if float(ret.rstrip())> 0: clashing_ligands.append(ligand_file_tfmd_path)
		ok_ligands = list(set(ligand_file_tfmd_paths)-set(clashing_ligands))
		if len(ok_ligands)==0: return None
		max_resnumber = find_last_res_number(compiled_ligand_file_path)
		if len(clashing_ligands): print "clashing:", clashing_ligands
	else:
		print "not exist?"
		ok_ligands    = ligand_file_tfmd_paths
		if len(ok_ligands)==0: return None
		max_resnumber = find_last_res_number(path+"/"+anchor)

	# clashing with the main chain?
	clashing_ligands = []
	for ligand_file_tfmd_path in ok_ligands:
		cmd = "{} {} {}".format(geom_overlap, ligand_file_tfmd_path, path+"/"+anchor)
		ret = subprocess.check_output(cmd, shell=True)
		if float(ret.rstrip())> 0: clashing_ligands.append(ligand_file_tfmd_path)

	remaining_ligands = sorted(list(set(ok_ligands) - set(clashing_ligands)))
	if len(remaining_ligands)==0: return None

	# compile what we have so far; renumber
	outfile = open (compiled_ligand_file_path,"a")
	new_ligands = []
	for ligand_file_tfmd_path in remaining_ligands:
		infile = open(ligand_file_tfmd_path,"r")
		prev_resnumber = -1
		output_resno = max_resnumber
		for line in infile:
			# -1 will grab the chain id too
			reslabel = line[res_number_pos-1:res_number_pos+res_number_length+1].replace(" ","")
			resnumber = int(reslabel[1:])
			if resnumber != prev_resnumber:
				output_resno += 1
				prev_resnumber = resnumber
				new_ligands.append(line[res_name_pos:res_name_pos+res_name_length].replace(" ",""))
			outfile.write(line[:res_number_pos] + "%4d"%output_resno + line[res_number_pos+4:])
		infile.close()
	outfile.close()

	new_ligands = list(set(new_ligands))

	return new_ligands

#########################################
def epitope2dist_string(pdb1, pdb2):
	cmd = "{} {} {} 10.0".format(geom_epitope, pdb1,pdb2)
	distances = {}
	for line in subprocess.check_output(cmd, shell=True).split("\n"):
		field = line.lstrip().rstrip().split()
		if len(field)<2: continue
		#chain = field[0][0] # not using that info
		resn  = field[0][1:4]
		distance = float(field[1])
		distances[resn] = distance
	return distances

#########################################
def split_into_compounds(compound_file):
	infile = open (compound_file, "r")
	outfile = {}
	for line in infile:
		if line[:6]!='HETATM': continue
		resname = line[res_name_pos:res_name_pos+res_name_length].replace(' ','')
		res_number = line[res_number_pos:res_number_pos+res_number_length].replace(' ','')
		chain = line[chain_pos]
		res_key = "_".join([resname,res_number,chain])
		if not outfile.has_key(res_key):
			outf_name= compound_file.replace('.pdb','')+".{}.pdb".format(res_key)
			outfile[res_key] = open(outf_name,"w")
		outfile[res_key].write(line)

	for outf in outfile.values():outf.close()
	infile.close()

	return [outf.name for outf in outfile.values()]

#########################################
def strip_and_glue (main_model_info, compiled_ligand_file_path, ligand_list, scratch):

	main_model_path, main_model, chains_in_main_model = main_model_info
	# the last field in the name is the source (pdb or swissmodel) which we'll replace with "compiled"
	compiled_model = None
	for label in ['peptide','swissmodel','ivana']:
		if label in main_model:
			compiled_model = main_model.replace(label,"compiled")
	if not compiled_model:
		print "name convention changed?"
		exit()
	os.chdir(scratch)

	outfile = open(compiled_model,"w")
	infile  = open(main_model_path+"/"+main_model, "r")
	for line in infile:
		if line[:4]=='ATOM' or line[:6] == 'HETATM':
			resname = line[res_name_pos:res_name_pos+res_name_length].replace(' ','')
			if resname in crystallographic_additives + ligand_list: continue
		outfile.write(line)
	infile.close()
	outfile.close()

	#concatenate protein and ligand files
	cmd = "cat {} >> {}".format(compiled_ligand_file_path, compiled_model)
	subprocess.call(cmd, shell=True)

	#extract main chain into one file, and the rest in the other
	# extract peptide only for chain chains_in_main_model[0]
	cmd = "{} {} -p -c{} >> {}".format(extract_chain,compiled_model,chains_in_main_model[0], "mainchain.pdb")
	subprocess.call(cmd, shell=True)

	distances = {} # the dictionary will contain distances to each of the remarkable points - interface to other chains, ligands
	for ligand_filename in split_into_compounds(compiled_ligand_file_path):
		compound_key = ligand_filename.split(".")[-2]
		distances[compound_key] = epitope2dist_string("mainchain.pdb",ligand_filename)

	for chain in chains_in_main_model[1:]:
		chain_pdb = "chain{}.pdb".format(chain)
		cmd = "{} {} -p -c{} >> {}".format(extract_chain,compiled_model,chain,chain_pdb)
		subprocess.call(cmd, shell=True)
		distances["chain"+chain] = epitope2dist_string("mainchain.pdb",chain_pdb)

	distance_string = json.dumps(distances)
	return compiled_model, distance_string

#########################################
def prepare_main_model(swissmodel_dir, model,scratch):
	cmd = "find {} -name {}".format(swissmodel_dir,model)
	path = subprocess.check_output(cmd, shell=True).rstrip()
	model_name = model.replace("swissmodel","peptide")
	modelfile = scratch+"/"+model_name
	inf = open(path,"r")
	outf = open(modelfile,"w")
	chains = set()
	for line in inf:
		if line[:6]=='HETATM' and line[chain_pos]=="_": continue # strip swissmodel ligands
		if line[:4]=='ATOM': chains.add(line[chain_pos])
		if line[:3]=='END': continue # we are going to glue the ligands at the bottom of the file
		outf.write(line)
	inf.close()
	outf.close()
	chains = sorted(list(chains))
	return scratch, model_name, chains

#########################################
def compile_model(cursor, gene_symbol,scratch):

	chains = {}
	distance_string = ""

	qry = "select  distinct(main_model) from model_elements where gene_symbol='%s'" % gene_symbol
	ret = search_db(cursor,qry)
	if len(ret)>2:
		print "different models for the same gene?"
		exit()
	main_model_info = prepare_main_model(swissmodel_dir, ret[0][0], scratch)
	main_model_path = main_model_info[0]
	main_model      = main_model_info[1]

	qry  = "select  pdb_chain, pdb_ligand, ligand_function, pdb_pct_identical_to_uniprot from model_elements"
	qry += " where gene_symbol='%s' order by pdb_pct_identical_to_uniprot desc"  % gene_symbol
	ret = search_db (cursor,qry)
	compiled_ligands_file_path = scratch + "/compiled_ligands.pdb"
	compiled_ligand_list = []
	ligand_functions = {}
	for line in ret:
		pdb_chain, pdb_ligand, ligand_function, pdb_pct_identical_to_uniprot = line
		if not ligand_functions.has_key(pdb_chain): ligand_functions[pdb_chain]= {}
		ligand_functions[pdb_chain][pdb_ligand] = ligand_function
		cmd = "{} {}".format(pdbdown, pdb_chain[:-1]) # this will check if it already exists
		subprocess.call(cmd, shell=True)
		ligand_file_path, new_chain = extract_ligand(pdb_path, scratch, pdb_chain[:-1]+".pdb", pdb_chain[-1], pdb_ligand)
		if os.path.getsize(ligand_file_path)==0: continue
		# use the transformation matrix to map all ligands to the main_model structure
		if new_chain:
			pdbchain = new_chain
		else:
			pdbchain = pdb_chain[-1]
		ligand_file_tfmd_paths = map_ligand_to_model (main_model_info,pdb_path,pdb_chain[:-1]+".pdb", pdbchain, ligand_file_path, scratch)
		if not ligand_file_tfmd_paths or len(ligand_file_tfmd_paths)==0: continue
		# remove clashing ligands adn ligands far from the main chain
		# add ligands which are not clashing with the existing ones
		new_ligands = merge_ligands(main_model_path, main_model, compiled_ligands_file_path, ligand_file_tfmd_paths, scratch)
		if new_ligands:
			compiled_ligand_list += new_ligands
			compiled_ligand_list  = list(set(compiled_ligand_list))

	print compiled_ligand_list


	if os.path.exists(compiled_ligands_file_path) and os.path.getsize(compiled_ligands_file_path)>0:
		chains = main_model_info[2]
		compiled_model, distance_string = strip_and_glue(main_model_info, compiled_ligands_file_path,
		                                                              compiled_ligand_list, scratch)
		# move out of struct_scratch
		if compiled_model and len(compiled_model)>0:
			os.rename(compiled_model,  main_model_path+"/"+compiled_model)
	else:
		compiled_model = None

	#print compiled_ligand_list
	#print distance_string
	return main_model_path, compiled_model, chains, compiled_ligand_list, distance_string, ligand_functions

##########################################
def store_ligands(cursor, approved_symbol, model_path, compiled_model, chains, compiled_ligands, distance_string):

	final_compiled_name = approved_symbol+"_"+compiled_model
	shutil.copy(model_path +"/" + compiled_model, compiled_model_repository + "/" + final_compiled_name)
	print "copied {}  to {}".format(model_path +"/" + compiled_model, compiled_model_repository + "/" + final_compiled_name)
	fixed_fields  = {'gene_symbol':approved_symbol, 'structure_file': final_compiled_name}
	chain = chains[0]
	other_chains = ",".join(chains[1:])
	# TODO: FE2 is an ion, and is 3 letters long ...
	ions = ",".join(list(set([x for x in compiled_ligands if x in physiological_ions])))
	substrates = ",".join(list(set([x for x in compiled_ligands if not x in physiological_ions])))
	update_fields = {'chain': chain, 'other_chains': other_chains, 'ions': ions,
	                 'substrates': substrates, 'distance_to_ligands':distance_string}
	store_or_update (cursor, 'monogenic_development.structures', fixed_fields, update_fields, verbose=False)

	return

##########################################
def main():

	for dependency in [swissmodel_dir, scratch_dir, struct, pdb_affine,
	                   geom_overlap, extract_chain, compiled_model_repository,
	                   pdb_path, pdbdown, pdb_chain_rename]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	db, cursor = connect()
	switch_to_db(cursor, "monogenic_development")
	qry = "select distinct(gene_symbol) from model_elements"
	ret = search_db(cursor,qry)

	for line in ret:
		gene_symbol = line[0]
		if gene_symbol!='PAH': continue
		print gene_symbol
		os.chdir(cwd)
		scratch = scratch_dir + "/" + gene_symbol
		shutil.rmtree(scratch, ignore_errors=True)
		os.makedirs(scratch)
		model_path, compiled_model, chains, compiled_ligands, distance_string, ligand_functions = compile_model(cursor,gene_symbol,scratch)
		if not compiled_model: continue
		print compiled_model, chains, compiled_ligands
		print distance_string
		# store compiled model in hte pdb directory of the monogenic server
		# store the list  of the ligands to the database
		store_ligands (cursor, gene_symbol, model_path, compiled_model, chains, compiled_ligands, distance_string)
		#shutil.rmtree(scratch, ignore_errors=True)



	########################################
if __name__ == '__main__':
	main()
