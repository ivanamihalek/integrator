
import urllib2
from bs4   import BeautifulSoup

##########################################
def ligands_from_pdb(pdb_id):
	pdb_request = "https://www.rcsb.org/pdb/rest/ligandInfo?structureId=%s" % pdb_id
	response = urllib2.urlopen(pdb_request)
	html = response.read()
	soup = BeautifulSoup(html, 'html.parser')

	return_list = []

	for ligand in soup.find_all('ligand'):
		ligand_info = {}
		ligand_info['chemical_id'] = str(ligand['chemicalid'])
		ligand_info['ligand_type'] = str(ligand['type'])
		ligand_info['chemical_name'] = str(ligand.chemicalname.string)
		ligand_info['smiles'] =  str(ligand.smiles.string)
		return_list.append(ligand_info)

	return return_list

##########################################
def uniprot_from_pdb_chain(pdb_chain_id):

	pdb_request = "https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=%s.%s" % (pdb_chain_id[:-1],pdb_chain_id[-1])
	response = urllib2.urlopen(pdb_request)
	html = response.read()
	soup = BeautifulSoup(html, 'html.parser')
	uniprot_id = str(soup.find(dbsource="UniProt")['dbaccessionid'])
	return uniprot_id

##########################################
def ec_cofactors_from_uniprot(uniprot_id):

	uniprot_request = "http://www.uniprot.org/uniprot/%s.xml" % uniprot_id
	response = urllib2.urlopen(uniprot_request)
	html = response.read()
	soup = BeautifulSoup(html, 'html.parser')
	ec_numbers = ";".join([str(ec['id']) for ec in soup.find_all('dbreference',type="EC")])
	cofactors = []
	for cof in soup.find_all('comment',type="cofactor"):
		for c in  cof.findChildren():
			if c.name=='name': str(cofactors.append(c.text))
	return ec_numbers, ";".join(cofactors)