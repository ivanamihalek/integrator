#!/usr/bin/python
from integrator_utils.mysql import *
import os

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
def pdb_to_dict_of_chars(path,filename):
	sequence = {}
	chain_pos =    21  #length 1
	res_name_pos = 17
	res_name_length = 3
	res_number_pos = 22
	res_number_length = 4
	aa_translation = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
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
	pdb_sequence    = pdb_to_dict_of_chars(path,model)
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
def extract_ligands(model):
	return
#########################################
def align (model, anchor):
	return
#########################################
def transform(ligand_file, tfm_matrix_file):
	return
#########################################
def compile_ligands():
	return

#########################################
def strip_and_glue (anchor, compiled_ligand_file):
	return

#########################################
def compile_model( model_info, path, protein):
	print "compiling model for", protein
	# find the anchor models
	# there can be multiple models if they are non-overlapping
	# two models are non-overlapping if they overlap by less than 80%
	# assign all other models to one of the anchors
	cluster = find_model_clusters(model_info, path)
	for anchor, anchored in cluster.iteritems():
		print anchor, anchored
	# aligns structurally all models to their anchor
	compiled_ligand_list = []
	for anchor, anchored in cluster.iteritems():
		# lignads do no have a name that signs like a residue or modified residue
		# ligand is not water and is not a crystallographic additive
		anchor_ligand_file = extract_ligands(model)
		for model in anchored:
			ligands, ligand_file = extract_ligands(model)
			compiled_ligand_list += ligands
			if not ligand_file: continue # no ligands here
			# use the transformation matrix to map all ligands to the anchor structure
			tfm_matrix_file = align (model, anchor)
			ligand_file_tfmd = transform(ligand_file, tfm_matrix_file)
			# if any ligand overlaps with one already present in the anchor, drop it
			# otherwise put it all in a new compiled ligand file
			compiled_ligand_file = compile_ligands()
		# strip anchor of all ligands
		# put in the compiled ligand files - call the whole thing compiled_from_to
		compiled_model = strip_and_glue (anchor, compiled_ligand_file)
	return compiled_model, compiled_ligand_list

##########################################
def store_ligands (cursor, compiled_ligand_list):
	return
##########################################
def distance_to_ligands(compiled_model, compiled_ligand_list):
	return
##########################################
def store_distance_string (cursor, distance_string):
	return

##########################################
def main():
	swissmodel_dir = "/databases/swissmodel"
	swissmodel_meta_file = swissmodel_dir+"/index_human.csv"
	model_info = parse_meta(swissmodel_meta_file) # there are more files (proteins/genes) here than qhat we will be using
	db, cursor = connect()


	#switch_to_db(cursor, "monogenic_development")
	qry = "select * from monogenic_development.diseases"
	ret = search_db(cursor, qry)
	for line in ret:
		[id, name_short, name_long, omim_ids, description, prim, sec] = line
		print "="*60
		print name_short
		for omim_id in omim_ids.split(";"):
			qry = "select approved_symbol from omim_genemaps where mim_number='%s'" % omim_id
			ret2 = search_db(cursor, qry)
			for line2 in ret2:
				protein = line2[0]
				print "\t", protein
				path = "/". join([swissmodel_dir,protein[0],protein])
				if not os.path.exists(path):
					print "\t", "\t", path, "not found"
					continue
				swissfound = False
				for model in next(os.walk(path))[2]:
					if 'swissmodel' in model:
						swissfound = True
						# the percentage of positions in each chain that is identical to the numberin the uniprot sequence
						identical_pct = check_res_numbers(cursor, swissmodel_dir, path, model)
						model_info[model]['identical_pct'] = identical_pct
				if not swissfound:
					print "\t", "\t", 'swissmodel not found'
					continue #  sometimes the models is not found - not sure if it is an old index file

				compiled_model, compiled_ligand_list = compile_model(model_info, path,protein)
				# store compiled model in hte pdb directory of the monogenic server
				# store the list  of the ligands to the database
				store_ligands (cursor, compiled_ligand_list)
				# measure the distance of all residues to ligands - store in the same order as lignads themselves
				distance_string = distance_to_ligands(compiled_model, compiled_ligand_list)
				store_distance_string (cursor, distance_string)

 ########################################
if __name__ == '__main__':
	main()
