#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess

swissmodel_dir = "/databases/swissmodel"
pdb_path       = "/databases/pdb/structures"
scratch_dir    = "/home/ivana/scratch"
struct         = "/home/ivana/projects/enzyme_modeling/code/struct/struct"
pdb_affine     = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_affine_tfm.pl"
pdbdown        = "/home/ivana/pypeworks/integrator/integrator_utils/pdbdownload.pl"
extract_chain  = "/home/ivana/pypeworks/integrator/integrator_utils/pdb_extract_chain.pl"
geom_epitope   = "/home/ivana/pypeworks/integrator/24_pdb/integrator_utils/geom_epitope.pl"
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
def extract_ligand(path, scratch, filename, chain, ligand_resn):

	outfilename = "/".join([scratch, filename.replace('.pdb','')+".{}.{}.pdb".format(chain,ligand_resn)])
	outfile = open(outfilename,"w")
	infile  = open (path+"/"+filename, "r")
	for line in infile:
		if line[:6] != 'HETATM': continue
		if line[chain_pos] != chain: continue
		resname = line[res_name_pos:res_name_pos+res_name_length].replace(' ','')
		if resname != ligand_resn: continue
		outfile.write(line)
	infile.close()
	outfile.close()
	return outfilename

#########################################
def map_ligand_to_model(main_model_info,pdb_path,ligand_containing_structure, ligand_chain, ligand_file_path, scratch):

	struct_scratch = scratch + "/" + "struct_scratch"
	shutil.rmtree(struct_scratch,ignore_errors=True)
	os.makedirs(struct_scratch)
	os.chdir(struct_scratch)

	# find the  first chain in anchor
	[path, main_model, main_model_chains] = main_model_info

	# do we have the tfm already by any chance
	tfm_key = "{}_{}_{}".format(main_model,main_model_chains[0],ligand_containing_structure)
	if transform.has_key(tfm_key):
		outf = open("tmp.tfm","w")
		outf.write(transform[tfm_key])
		outf.close()
	else:
		from_structure_path = pdb_path+"/"+ligand_containing_structure
		to_structure_path   = path+"/"+main_model
		cmd = "%s -from %s -c1 %s -to %s -c2 %s" % (struct, from_structure_path, ligand_chain, to_structure_path, main_model_chains[0])
		subprocess.call(cmd, shell=True)
		almt_file = "{}{}_to_{}{}.0.aln".format(ligand_containing_structure.replace('.pdb',''), ligand_chain, main_model.replace('.pdb',''),main_model_chains[0])
		if not os.path.exists(almt_file) or os.stat(almt_file).st_size == 0: return None
		print "blah 1"
		cmd = "grep tfm {} -A 3 | tail -n 3 | sed 's/{}//g' > tmp.tfm".format(almt_file,"%")
		subprocess.call(cmd, shell=True)
		if not os.path.exists('tmp.tfm') or os.stat('tmp.tfm').st_size == 0: return None
		print "blah 2"
		inf = open("tmp.tfm","r")
		transform[tfm_key] = inf.read()
		inf.close()

	ligand_file_tfmd_path = ligand_file_path.replace('.pdb','') + ".tfmd.pdb"
	cmd = "{} {} tmp.tfm > {}".format(pdb_affine, ligand_file_path, ligand_file_tfmd_path)
	subprocess.call(cmd, shell=True)

	os.chdir(scratch)
	shutil.rmtree(struct_scratch,ignore_errors=True)

	return ligand_file_tfmd_path

#########################################
def merge_ligands(path, anchor, compiled_ligand_file_path, ligand_file_tfmd_path, scratch):

	if not os.path.exists(compiled_ligand_file_path):
		shutil.copy(ligand_file_tfmd_path,compiled_ligand_file_path)
		return None

	os.chdir(scratch)
	# do I have anything closer than1 1 A?
	#note that I am taking the footprint on the ligand, to remove clashing molecules
	clashing_molecules = []
	cmd = "{} {} {} 1.0".format(geom_epitope, ligand_file_tfmd_path, path+"/"+anchor)
	ret = subprocess.check_output(cmd, shell=True)
	for line in  ret.split("\n"):
		if line.replace(" ","")=="": continue
		clashing_molecules.append(line[:4].replace(" ",""))

	cmd = "{} {} {} 1.0".format(geom_epitope, ligand_file_tfmd_path,compiled_ligand_file_path)
	ret = subprocess.check_output(cmd, shell=True)
	for line in  ret.split("\n"):
		if line.replace(" ","")=="": continue
		clashing_molecules.append(line[:4].replace(" ",""))

	# get rid of duplicates
	clashing_molecules = list(set(clashing_molecules))
	# remove clashing molecules from the file; renumber
	outfile = open ("tmp_ligand.pdb","w")

	# write the old file, find the last res number
	infile = open(compiled_ligand_file_path,"r")
	max_resnumber = -1
	for line in infile:
		resnumber = int(line[res_number_pos:res_number_pos+res_number_length])
		outfile.write(line)
		if max_resnumber<resnumber: max_resnumber=resnumber
	infile.close()

	infile = open(ligand_file_tfmd_path,"r")
	prev_resnumber = -1
	output_resno = max_resnumber
	new_ligands = []
	for line in infile:
		# -1 will grab the chain id too
		reslabel = line[res_number_pos-1:res_number_pos+res_number_length+1].replace(" ","")
		if reslabel in clashing_molecules: continue
		resnumber = int(reslabel[1:])
		if resnumber != prev_resnumber:
			output_resno += 1
			prev_resnumber = resnumber
			new_ligands.append(line[res_name_pos:res_name_pos+res_name_length].replace(" ",""))
		outfile.write(line[:res_number_pos] + "%4d"%output_resno + line[res_number_pos+4:])
	infile.close()
	outfile.close()
	os.rename("tmp_ligand.pdb",compiled_ligand_file_path)
	return new_ligands

#########################################
def strip_and_glue (main_model_info, compiled_ligand_file_path, ligand_list, scratch):

	path, anchor, chains_in_anchor = main_model_info
	# the last field in the name is the source (pdb or swissmodel) which we'll replace with "compiled"
	compiled_model = "_".join(anchor.split("_")[:-1] + ["compiled.pdb"])
	os.chdir(scratch)

	outfile = open(compiled_model,"w")
	infile  = open(path+"/"+anchor, "r")
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
	# extract peptide only for chain chains_in_anchor[0]
	cmd = "{} {} -p -c{} >> {}".format(extract_chain,compiled_model,chains_in_anchor[0], "mainchain.pdb")
	subprocess.call(cmd, shell=True)

	# extract everything  except mainchain (inverse; -i) for chain chains_in_anchor[0]
	cmd = "{} {} -i -c{} >> {}".format(extract_chain,compiled_model,chains_in_anchor[0], "else.pdb")
	subprocess.call(cmd, shell=True)

	distances = [] # residues to any of the ligands
	cmd = "{} {} {} 10.0".format(geom_epitope, "mainchain.pdb","else.pdb")
	#cmd = "{} {} {} 10.0".format(geom_epitope, "mainchain.pdb",compiled_ligand_file_path)
	for line in subprocess.check_output(cmd, shell=True).split("\n"):
		field = line.lstrip().rstrip().split()
		if len(field)<2: continue
		chain = field[0][0]
		if chains_in_anchor[0]==chain:
			resn  = field[0][1:4]
			distance = float(field[1])
			distances.append("%s:%.2f"%(resn,distance))
	distance_string = ",".join(distances)

	return compiled_model, distance_string

#########################################
def find_main_model_info(swissmodel_dir,model):
	cmd = "find {} -name {}".format(swissmodel_dir,model)
	path = subprocess.check_output(cmd, shell=True).rstrip()
	inf = open(path,"r")
	chains = set()
	for line in inf:
		if line[:4]=='ATOM': chains.add(line[chain_pos])
	inf.close()
	chains = sorted(list(chains))
	aux    = path.split("/")
	aux.pop()
	path = "/".join(aux)
	return path, model, chains

#########################################
def compile_model(cursor, gene_symbol):

	os.chdir(cwd)
	scratch = scratch_dir + "/" + gene_symbol
	shutil.rmtree(scratch,ignore_errors=True)
	os.makedirs(scratch)
	prev_main_model = None

	compiled_ligands_for_model = {}
	compiled_model = ""
	chains = {}
	distance_strings = {}

	qry = "select  distinct(main_model) from model_elements where gene_symbol='%s'" % gene_symbol
	ret = search_db(cursor,qry)
	if len(ret)>2:
		print "different models for the same gene?"
		exit()
	main_model = ret[0][0]
	main_model_info = find_main_model_info(swissmodel_dir,main_model)
	main_model_path = main_model_info[0]
	print main_model_info

	qry = "select  id, metacyc_ligand, ligand_tanimoto from model_elements where gene_symbol='%s'" % gene_symbol
	max_tanimoto = {}
	max_id       = {}
	for line in  search_db(cursor,qry):
		[id, metacyc_ligand, ligand_tanimoto] = line
		if not max_tanimoto.has_key(metacyc_ligand) or ligand_tanimoto > max_tanimoto[metacyc_ligand]:
			max_tanimoto[metacyc_ligand] = ligand_tanimoto
			max_id[metacyc_ligand]       = id


	compiled_ligands_file_path = scratch + "/compiled_ligands.pdb"
	compiled_ligand_list = []
	for metacyc_ligand, id_with_the_highest_tanimoto_score in max_id.iteritems():
		qry = "select  pdb_chain, pdb_ligand from model_elements where id=%d"  % id_with_the_highest_tanimoto_score
		ret  = search_db(cursor,qry)
		pdb_chain, pdb_ligand  = ret[0]
		cmd = "{} {}".format(pdbdown, pdb_chain[:-1]) # this will check if it already exists
		subprocess.call(cmd, shell=True)
		ligand_file_path = extract_ligand(pdb_path, scratch, pdb_chain[:-1]+".pdb", pdb_chain[-1], pdb_ligand)
		if os.path.getsize(ligand_file_path)==0: continue
		# use the transformation matrix to map all ligands to the main_model structure
		ligand_file_tfmd_path = map_ligand_to_model (main_model_info,pdb_path,pdb_chain[:-1]+".pdb", pdb_chain[-1], ligand_file_path, scratch)
		if not ligand_file_tfmd_path: continue
		# remove clashing ligands adn ligands far from the main chain
		# add ligands which are not clashing with the existing ones
		new_ligands = merge_ligands(main_model_path, main_model, compiled_ligands_file_path, ligand_file_tfmd_path, scratch)
		if not new_ligands: new_ligands=[pdb_ligand]
		compiled_ligand_list += new_ligands
		compiled_ligand_list  = list(set(compiled_ligand_list))

	if  os.path.exists(compiled_ligands_file_path) and  os.path.getsize(compiled_ligands_file_path)>0:
		# strip main_model of all ligands
		# put in the compiled ligand files - call the whole thing compiled_from_to
		# move out of scratch
		chains = main_model_info[2]
		compiled_model, distance_strings[main_model] = strip_and_glue(main_model_info, compiled_ligands_file_path, compiled_ligand_list, scratch)
		print compiled_model
		if compiled_model and len(compiled_model)>0:
			os.rename(compiled_model,  main_model_path+"/"+compiled_model)
	else:
		compiled_model = None

	shutil.rmtree(scratch,ignore_errors=True)

	distance_string =  ",".join(distance_strings.values())
	#print compiled_ligand_list
	#print distance_string
	return main_model_path, compiled_model, chains, compiled_ligand_list, distance_string

##########################################
def store_ligands(cursor, approved_symbol, path, compiled_model, chains, compiled_ligands, distance_string):

	shutil.copy(path +"/" + compiled_model, compiled_model_repository + "/" + compiled_model)
	fixed_fields  = {'gene_symbol':approved_symbol, 'structure_file': compiled_model}
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
	                   geom_epitope, extract_chain, compiled_model_repository,
	                   pdb_path, pdbdown]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	swissmodel_meta_file = swissmodel_dir + "/index_human.csv"

	db, cursor = connect()
	switch_to_db(cursor, "monogenic_development")
	qry = "select distinct(gene_symbol) from model_elements"
	ret = search_db(cursor,qry)
	for line in ret:
		gene_symbol = line[0]
		if gene_symbol!='MUT': continue
		print gene_symbol
		path, compiled_model, chains, compiled_ligands, distance_string = compile_model(cursor,gene_symbol)
		if not compiled_model: continue
		print compiled_model, chains, compiled_ligands
		print distance_string
		# store compiled model in hte pdb directory of the monogenic server
		# store the list  of the ligands to the database
		store_ligands (cursor, gene_symbol, path, compiled_model, chains, compiled_ligands, distance_string)

 ########################################
if __name__ == '__main__':
	main()
