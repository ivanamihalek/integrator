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
geom_overlap      = "/home/ivana/pypeworks/integrator/24_pdb/integrator_utils/geom_overlap.pl"
geom_epitope      = "/home/ivana/pypeworks/integrator/24_pdb/integrator_utils/geom_epitope.pl"
compiled_model_repository = "/home/ivana/monogenic/public/pdb"
cwd    = os.getcwd()

#########################################
chain_pos      = 21  #length 1
res_name_pos   = 17
res_name_length = 3
res_number_pos = 22
res_number_length = 4
aa_translation = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                  'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
#########################################
transform = {}
crystallographic_additives = ['HOH', "SO4", "GOL","PGO","PGR","EDO","EOH","DIO","SOG","HTG","CL","PEG","NAG", "SAM", "HEM","PE4","ACT","NA"]

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
		print tfm_key
		if transform.has_key(tfm_key):
			print transform[tfm_key]
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
	#shutil.rmtree(struct_scratch,ignore_errors=True)

	return ligand_file_tfmd_paths


##########################################

##########################################
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

	print "scratch", scratch
	print "main model", main_model_path+"/"+main_model
	print "compiled", compiled_model

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



########################################
def main():

	for dependency in [swissmodel_dir, scratch_dir, struct, pdb_affine,
	                   geom_overlap, extract_chain, compiled_model_repository,
	                   pdb_path, pdbdown, pdb_chain_rename]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	scratch = scratch_dir + "/PAH"
	ligand_containing_structure = "5fii.pdb"
	ligand_chain = "A"
	main_model_info = [scratch,"P00439_20_450_5fgj.1.A_ivana.pdb", ['A','B','C','D']]

	if False:
		ligand_path, something =  extract_ligand(pdb_path, scratch, ligand_containing_structure, ligand_chain, 'PHE')
		print ligand_path

		# [path, main_model, main_model_chains]
		map_ligand_to_model(main_model_info,pdb_path,ligand_containing_structure, ligand_chain, ligand_path, scratch)
	    # map_ligand_to_model(main_model_info,pdb_path,ligand_containing_structure, ligand_chain, ligand_file_path, scratch):

	compiled_ligand_list = ['H4B', 'DAH', 'PHE', 'LNR', 'FE2', 'FE', 'TIH', 'TRS']
	chains = main_model_info[2]
	compiled_ligands_file_path = scratch + "/compiled_ligands.pdb"
	compiled_model, distance_strings = strip_and_glue(main_model_info, compiled_ligands_file_path,
		                                                              compiled_ligand_list, scratch)
	print compiled_model
	print distance_strings

########################################
if __name__ == '__main__':
	main()
