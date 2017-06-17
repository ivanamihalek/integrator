#!/usr/bin/python
from integrator_utils.mysql import *
import os, shutil, subprocess

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
crystallographic_additives = ['HOH', "SO4", "GOL","PGO","PGR","EDO","EOH","DIO","SOG","HTG","CL","PEG","NAG","HEM"]
physiological_ions =["FE","FE2","MN","ZN","ZN2","MG","CU","CO","CD","MO","VA","NI","W", "SE","CA"]

##########################################
def parse_url(ur):
	# url of the form
	#https://swissmodel.expasy.org/repository/uniprot/P22304.pdb?from=34&to=550&template=5fql.1.A&provider=swissmodel
	[addr, fnm] = ur.split('?')
	uniprot_id = addr.split("/")[-1].replace('.pdb','')
	pdbname =  uniprot_id
	for field in fnm.split("&"):
		[k,v] = field.split("=")
		pdbname += "_" + v
	pdbname += ".pdb"
	return pdbname

##########################################
def pdb_to_sequence(path, filename):
	sequence = {}
	#open file
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

########################################
def check_res_numbers(cursor, swissmodel_dir, path, model):

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

########################################
def fix_by_hand(url):
	# some url's are broken
	manual_hack = {}
	#manual_hack['https://swissmodel.expasy.org/repository/uniprot/P35914.pdb?from=28&to=323&template=3mp3.1.A&provider=swissmodel'] =\
	#          "P35914_28_189_3mp3.1.B_swissmodel.pdb"
	if not manual_hack.has_key(url): return None, None
	fnm = manual_hack[url]
	[uni, start, end, template, provider] = fnm.replace('.pdb','').split("_")
	return fnm, {'start':start, 'end':end, 'qmean':0.0,'qmean_norm':0.0, 'provider':provider}

########################################
def parse_meta(swissmetafile):
	infile = open(swissmetafile,"r")
	model_info = {}
	for line in infile:
		if line[0]=='#': continue
		line = line.rstrip()
		[uniprot_id, coordinate_id, provider, start, end, template, qmean, qmean_norm,url] = line.split('\t')
		pdbname, hack = fix_by_hand(url)
		if pdbname:
			model_info[pdbname] = hack
		else:
			pdbname = parse_url(url)
			model_info[pdbname] = {'start':start, 'end':end, 'qmean':qmean,'qmean_norm':qmean_norm, 'provider':provider}
	return model_info

########################################
def overlap(model_info, A, B):
	
	As = int(model_info[A]['start'])
	Ae = int(model_info[A]['end'])
	Bs = int(model_info[B]['start'])
	Be = int(model_info[B]['end'])	
	
	if Ae<Bs or  Be<As:
		return None, None # A and B do not cover each other
	elif  As<=Bs  and  Be<=Ae:
		return A, B # A covers B
	elif  Bs<=As  and  Ae<=Be:
		return B, A # B covers A
	# from this point on we improvise
	# say, if difference is < 20 res, the longer one covers; otherwis, use both
	lenA = Ae - As
	lenB = Be - Bs
	len_overlap = float(min(Ae,Be) - max(As,Bs))
	#print A, B
	#print lenA, lenB, len_overlap
	if lenA>=lenB:
		if len_overlap/lenB>0.8:
			return A, B # A covers B
	elif len_overlap/lenA>0.8:
			return B, A
	return None, None

#########################################
def find_model_clusters(model_info, path):
	cluster = {}
	not_eligible_for_anchor = []
	for model in next(os.walk(path))[2]:
		# ignore junk in the directory
		if "pdb" != model[-3:] or 'compiled' in model: continue
		if model_info[model]['provider']!='swissmodel':
			not_eligible_for_anchor.append(model)
		else:
			if  model_info[model]['identical_pct']['A'] < 80:
				not_eligible_for_anchor.append(model)
				continue
			covering = None
			for anchor in cluster.keys():
				covering, covered = overlap(model_info, model, anchor)
				if covering: break
			if not covering:
				cluster[model] = []
			elif covered==model:
				cluster[covering].append(model)
			else:
				cluster[model] = cluster[covered][:]
				cluster[model].append(covered)
				del cluster[covered]
	for model in not_eligible_for_anchor:
		for anchor in cluster.keys():
			covering, covered = overlap(model_info, model, anchor)
			if covering: # if there is any overlpa, I'll take my chances
				cluster[anchor].append(model)
				break
	return cluster

#########################################
def extract_ligands(path, filename, scratch):

	outfilename = "/".join([scratch, filename.replace('.pdb','')+".ligands.pdb"])

	outfile = open(outfilename,"w")
	infile  = open (path+"/"+filename, "r")
	ligands = []
	for line in infile:
		if line[:6] != 'HETATM': continue
		chain = line[chain_pos]
		# I will take a peal of faith here and trut that Torsten and his buddies have maeked the
		# the ligands (as opposed to modified residues) by setting the chain to "_"
		if chain != "_": continue
		resname = line[res_name_pos:res_name_pos+res_name_length].replace(' ','')
		if resname in crystallographic_additives: continue
		if not resname in ligands: ligands.append(resname)
		outfile.write(line)
	infile.close()
	outfile.close()

	return ligands, outfilename

#########################################
def map_ligands_to_model (model_info, path, model, anchor, ligand_file_path, scratch):
	struct_scratch = scratch + "/" + "struct_scratch"
	shutil.rmtree(struct_scratch,ignore_errors=True)
	os.makedirs(struct_scratch)
	os.chdir(struct_scratch)

	# find the alphabetically first chain in anchor
	first_chain = sorted(model_info[anchor]['identical_pct'].keys())[0]
	cmd = "%s -from %s -to %s -c2 %s" % (struct, path+"/"+model, path+"/"+anchor, first_chain)
	subprocess.call(cmd, shell=True)
	almt_file = "{}_to_{}{}.0.aln".format(model.replace('.pdb',''),anchor.replace('.pdb',''),first_chain)
	if not os.path.exists(almt_file) or os.stat(almt_file).st_size == 0: return None
	cmd = "grep tfm {} -A 3 | tail -n 3 | sed 's/{}//g' > tmp.tfm".format(almt_file,"%")
	subprocess.call(cmd, shell=True)
	if not os.path.exists('tmp.tfm') or os.stat('tmp.tfm').st_size == 0: return None
	# TDOD check that the interesting ligands are really mapping to chain A, and not some place else
	ligand_file_tfmd_path = ligand_file_path.replace('.pdb', '') + ".tfmd.pdb"

	cmd = "{} {} tmp.tfm > {}".format(pdb_affine, ligand_file_path, ligand_file_tfmd_path)
	subprocess.call(cmd, shell=True)

	os.chdir(scratch)
	shutil.rmtree(struct_scratch,ignore_errors=True)

	return ligand_file_tfmd_path

#########################################
def merge_ligands(path, anchor, compiled_ligand_file_path, ligand_file_tfmd_path, scratch):
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
def strip_and_glue (path, anchor, chains_in_anchor, compiled_ligand_file_path, ligand_list, scratch):
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
def compile_model (model_info, path, protein):
	# find the anchor models
	# there can be multiple models if they are non-overlapping
	# two models are non-overlapping if they overlap by less than 80%
	# assign all other models to one of the anchors
	cluster = find_model_clusters(model_info, path)
	# aligns structurally all models to their anchor
	compiled_ligands_for_model = {}
	compiled_models = []
	chains = {}
	distance_strings = {}
	for anchor, anchored in cluster.iteritems():
		print anchor, anchored
		os.chdir(cwd)
		scratch = scratch_dir + "/" + anchor.replace('.pdb','')
		shutil.rmtree(scratch,ignore_errors=True)
		os.makedirs(scratch)
		# ligands do not have a name that sounds like a residue or modified residue
		# ligand is not water and is not a crystallographic additive
		anchor_ligands, anchor_ligand_file_path = extract_ligands(path, anchor, scratch)
		compiled_ligand_list      = anchor_ligands
		compiled_ligand_file_path = scratch + "/compiled_ligands.pdb"
		shutil.copy(anchor_ligand_file_path,compiled_ligand_file_path)
		for model in anchored:
			ligands, ligand_file_path = extract_ligands(path, model,scratch)
			if not ligands or len(ligands)==0: continue # no ligands here
			# use the transformation matrix to map all ligands to the anchor structure
			ligand_file_tfmd_path = map_ligands_to_model (model_info, path, model, anchor, ligand_file_path, scratch)
			if not ligand_file_tfmd_path: continue
			# remove clashing ligands adn ligands far from the main chain
			# add ligands which are not clashing with the existing ones
			new_ligands = merge_ligands(path, anchor, compiled_ligand_file_path, ligand_file_tfmd_path, scratch)
			compiled_ligand_list += new_ligands
		compiled_ligand_list = list(set(compiled_ligand_list))
		# strip anchor of all ligands
		# put in the compiled ligand files - call the whole thing compiled_from_to
		# move out of scratch
		chns = sorted(model_info[anchor]['identical_pct'].keys())
		compiled_model, distance_strings[anchor] = strip_and_glue (path, anchor, chns, compiled_ligand_file_path, compiled_ligand_list, scratch)
		chains[compiled_model] = chns
		if compiled_model and len(compiled_model)>0:
			compiled_models.append(compiled_model)
			compiled_ligands_for_model[compiled_model] = compiled_ligand_list
			os.rename(compiled_model,  path+"/"+compiled_model)

		shutil.rmtree(scratch,ignore_errors=True)

	distance_string =  ",".join(distance_strings.values())
	print compiled_ligand_list
	print distance_string
	return compiled_models, chains, compiled_ligands_for_model, distance_string

##########################################
def store_ligands (cursor, approved_symbol, path, compiled_models, chains, compiled_ligands_for_model, distance_string):

	for model in compiled_models:
		shutil.copy(path+"/"+model, compiled_model_repository+"/"+model)
		fixed_fields  = {'gene_symbol':approved_symbol, 'structure_file': model}
		chain = chains[model][0]
		other_chains = ",".join(chains[model][1:])
		# TODO: FE2 is an ion, and is 3 letters long ...
		ions = ",".join(list(set([x for x in compiled_ligands_for_model[model] if x in physiological_ions] )))
		substrates = ",".join(list(set([x for x in compiled_ligands_for_model[model] if not x in physiological_ions])))
		update_fields = {'chain': chain, 'other_chains': other_chains, 'ions': ions,
		                 'substrates': substrates, 'distance_to_ligands':distance_string}
		store_or_update (cursor, 'monogenic_development.structures', fixed_fields, update_fields,verbose=False)

	return

##########################################
def main():

	for dependency in [swissmodel_dir, scratch_dir, struct, pdb_affine,
	                   geom_epitope, extract_chain, compiled_model_repository]:
		if os.path.exists(dependency): continue
		print dependency, " not found"
		exit()

	swissmodel_meta_file = swissmodel_dir+"/index_human.csv"
	model_info = parse_meta(swissmodel_meta_file) # there are more files (proteins/genes) here than what we will be using

	db, cursor = connect()
	#switch_to_db(cursor, "monogenic_development")
	qry = "select * from monogenic_development.diseases"
	ret = search_db(cursor, qry)
	for line in ret:
		[id, name_short, name_long, omim_ids, description, prim, sec] = line
		#print "="*60
		#print name_short
		for omim_id in omim_ids.split(";"):
			qry = "select approved_symbol from omim_genemaps where mim_number='%s'" % omim_id
			ret2 = search_db(cursor, qry)
			for line2 in ret2:
				approved_symbol = line2[0]
				if approved_symbol != 'CBS': continue
				print "\t", approved_symbol
				path = "/". join([swissmodel_dir,approved_symbol[0],approved_symbol])
				if not os.path.exists(path):
					print "\t", "\t", path, "not found"
					continue
				swissfound = False
				for model in next(os.walk(path))[2]:
					if 'compiled' in model: continue
					if 'swissmodel' in model: swissfound = True
					# the percentage of positions in each chain that is identical to the number in the uniprot sequence
					identical_pct = check_res_numbers(cursor, swissmodel_dir, path, model)
					model_info[model]['identical_pct'] = identical_pct
				if not swissfound:
					print "\t", "\t", 'swissmodel not found'
					continue #  sometimes the models is not found - not sure if it is an old index file

				compiled_models, chains, compiled_ligands_for_model, distance_string = compile_model(model_info, path, approved_symbol)
				#print compiled_models, chains, compiled_ligands_for_model
				#print distance_string
				# store compiled model in hte pdb directory of the monogenic server
				# store the list  of the ligands to the database
				store_ligands (cursor, approved_symbol, path, compiled_models, chains, compiled_ligands_for_model, distance_string)

 ########################################
if __name__ == '__main__':
	main()
