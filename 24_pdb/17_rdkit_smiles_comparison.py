#!/usr/bin/env python
# the only way I coyuld istall rdkit was with conda, so run with
# source /home/ivana/miniconda2/bin/activate rdkit-env; ./17_rdkit_smiles_comparison.py ; source  /home/ivana/miniconda2/bin/deactivate rdkit-env;

# (I could not move thow whole thing to conda either, because I could not get it to talk with mysql

import sys
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


##########################################
def tanimoto (smiles1, smiles2):
	ms = [Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)]
	fps = [FingerprintMols.FingerprintMol(x) for x in ms]
	similarity = DataStructs.FingerprintSimilarity(fps[0],fps[1])
	return similarity

##########################################
def main():
	if len(sys.argv) < 3:
		print  "usage: %s  <smiles1> <smiles2>" % sys.argv[0]
	smiles1 = sys.argv[1]
	smiles2 = sys.argv[2]
	print tanimoto (smiles1, smiles2)
	return

#########################################
if __name__ == '__main__':
	main()
