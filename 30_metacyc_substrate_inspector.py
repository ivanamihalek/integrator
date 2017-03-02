#!/usr/bin/python

# the core table in metacyc seems to be enzrxns
# it uses field | to link to

# enzyme                   | protein     --->  components | protein
# reaction                 | reaction    --->  rxn_left and rxn_right | compound
# cofactors                | compound
# required_protein_complex | protein (?)
# regulated_by             | regulation
# alternative_substrates   | compound
# alternative_cofactors    | compound


from integrator_utils.mysql import *
import os

##########################################
def get_enzrxn_ids(db, cursor):
	qry  = "select from_id from  metacyc_edges where from_field='enzrxn'"
	rows = search_db(cursor, qry)
	if not rows or len(rows)==0:
		print "no human enrxns found (?!)"
		cursor.close()
		db.close()
		exit()

	return [row[0] for row in rows]

##########################################
def get_substrates (cursor, reaction_id):
	subs = []
	qry  = "select rxn_left, rxn_right from metacyc_reactions where unique_id like '%s'" % reaction_id
	rows = search_db(cursor, qry)
	if not rows or len(rows)==0: return subs
	for row in rows:
		for entry in row:
			subs += entry.replace(' ','').split(";")

	return subs

##########################################
def main():
	db, cursor = connect()

	reaction_substrates = []
	enzrxn_ids = get_enzrxn_ids(db, cursor)
	for enzrxn_id in enzrxn_ids:
		# also process alternative substrates here
		qry  = "select reaction from metacyc_enzrxns where unique_id='%s'" % enzrxn_id
		rows = search_db(cursor, qry)
		if not rows or len(rows)==0:
			print "no reaction for %d ?!" % enzrxn_id
			continue
		for row in rows:
			reaction_substrates +=  get_substrates (cursor,  row[0])
	histogram = {}
	for sub in set(reaction_substrates): histogram[sub] = 0
	for sub in reaction_substrates: histogram[sub] += 1
	sorted_hist = sorted(histogram.items(), key=lambda x:x[1]) #list of tuples ordered by value
	print len(sorted_hist)
	for subs, occurence in sorted_hist:
		print "  %4d   %s  "  % (occurence, subs)

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()


############################################
# further problems
# the name of the compound seems to be cut:
# mysql> select rxn_left from metacyc_reactions where rxn_left like '%PHOSPHATIDYL-MYO-INOSITOL%';
# +------------------------------------------------+
# | rxn_left                                       |
# +------------------------------------------------+
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; WATER |
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; WATER |
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; ATP   |
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; WATER |
# | GDP-MANNOSE; 1-PHOSPHATIDYL-MYO-INOSITOL       |
# +------------------------------------------------+
# 5 rows in set (0.01 sec)
#

# but then id does not actually appear in the compound table ....

# mysql> select unique_id  from metacyc_compounds where unique_id like '%PHOSPHATIDYL-MYO-INOSITOL%';
# +---------------------------------------+
# | unique_id                             |
# +---------------------------------------+
# | 1-PHOSPHATIDYL-MYO-INOSITOL-2-MANNOSE |
# | 1-PHOSPHATIDYL-MYO-INOSITOL           |
# +---------------------------------------+
# 2 rows in set (0.01 sec)
#
# mysql> select rxn_right from metacyc_reactions where rxn_right like '%PHOSPHATIDYL-MYO-INOSITOL%';
# +------------------------------------------------------+
# | rxn_right                                            |
# +------------------------------------------------------+
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; ADP; PROTON |
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; ADP; PROTON |
# | PROTON; 1-PHOSPHATIDYL-MYO-INOSITOL-2-MANNOSE; GDP   |
# | PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA; |Pi|        |
#+------------------------------------------------------+

# how common is this?

############################################
# the most common inout/output of a reaction
# there is a jump between 22 and 26 - take this as the cutoff
#     21   |Peptides-holder|
#     22   GTP
#     22   GLUTATHIONE
#     26   MALONYL-COA
#     29   PHOSPHATIDYL-MYO-INOSITOL-45-BISPHOSPHA
#     30   ADENOSYL-HOMO-CYS
#     31   S-ADENOSYLMETHIONINE
#     32   |Acceptor|
#     32   |Donor-H2|
#     32   |Red-NADPH-Hemoprotein-Reductases|
#     32   |Ox-NADPH-Hemoprotein-Reductases|
#     33   UDP-N-ACETYL-D-GLUCOSAMINE
#     38   HYDROGEN-PEROXIDE
#     38   3-5-ADP
#     39   PAPS
#     41   AMMONIUM
#     41   2-KETOGLUTARATE
#     43   NADH-P-OR-NOP
#     43   NAD-P-OR-NOP
#     51   GLT
#     52   AMP
#     82   PPI
#     85   ACETYL-COA
#     86   UDP
#    112   NADPH
#    113   CARBON-DIOXIDE
#    115   NADP
#    127   NADH
#    133   NAD
#    182   OXYGEN-MOLECULE
#    205   CO-A
#    212   ADP
#    235   |Pi|
#    259   ATP
#    666   WATER
#   1003   PROTON
