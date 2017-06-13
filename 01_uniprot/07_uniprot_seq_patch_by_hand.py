#! /usr/bin/python

from integrator_utils.mysql import *

sequences = {}
sequences['Q96NT5'] = {}
sequences['Q96NT5']['sequence'] = "MEGSASPPEKPRARPAAAVLCRGPVEPLVFLANFALVLQGPLTTQYLWHRFSADLGYNGT" \
                      "RQRGGCSNRSADPTMQEVETLTSHWTLYMNVGGFLVGLFSSTLLGAWSDSVGRRPLLVLA" \
                      "SLGLLLQALVSVFVVQLQLHVGYFVLGRILCALLGDFGGLLAASFASVADVSSSRSRTFR" \
                      "MALLEASIGVAGMLASLLGGHWLRAQGYANPFWLALALLIAMTLYAAFCFGETLKEPKST" \
                      "RLFTFRHHRSIVQLYVAPAPEKSRKHLALYSLAIFVVITVHFGAQDILTLYELSTPLCWD" \
                      "SKLIGYGSAAQHLPYLTSLLALKLLQYCLADAWVAEIGLAFNILGMVVFAFATITPLMFT" \
                      "GYGLLFLSLVITPVIRAKLSKLVRETEQGALFSAVACVNSLAMLTASGIFNSLYPATLNF" \
                      "MKGFPFLLGAGLLLIPAVLIGMLEKADPHLEFQQFPQSP"

sequences['Q96NT5']['ensembl_transcript_id'] = 'ENST00000612814'
sequences['Q96NT5']['ensembl_protein_id']    = 'ENSP00000480703'


sequences['O43826'] = {}
sequences['O43826']['sequence']  ="MAAQGYGYYRTVIFSAMFGGYSLYYFNRKTFSFVMPSLVEEIPLDKDDLGFITSSQSAAY" \
                      "AISKFVSGVLSDQMSARWLFSSGLLLVGLVNIFFAWSSTVPVFAALWFLNGLAQGLGWPP" \
                      "CGKVLRKWFEPSQFGTWWAILSTSMNLAGGLGPILATILAQSYSWRSTLALSGALCVVVS" \
                      "FLCLLLIHNEPADVGLRNLDPMPSEGKKGSLKEESTLQELLLSPYLWVLSTGYLVVFGVK" \
                      "TCCTDWGQFFLIQEKGQSALVGSSYMSALEVGGLVGSIAAGYLSDRAMAKAGLSNYGNPR" \
                      "HGLLLFMMAGMTVSMYLFRVTVTSDSPKLWILVLGAVFGFSSYGPIALFGVIANESAPPN" \
                      "LCGTSHAIVGLMANVGGFLAGLPFSTIAKHYSWSTAFWVAEVICAASTAAFFLLRNIRTK" \
                      "MGRVSKKAE"

sequences['O43826']['ensembl_transcript_id'] = 'ENST00000545985'
sequences['O43826']['ensembl_protein_id']    = 'ENSP00000475241'


sequences['O14841'] = {}
sequences['O14841']['sequence'] = "MGSPEGRFHFAIDRGGTFTDVFAQCPGGHVRVLKLLSEDPANYADAPTEGIRRILEQEAG" \
                      "MLLPRDQPLDSSHIASIRMGTTVATNALLERKGERVALLVTRGFRDLLHIGTQARGDLFD" \
                      "LAVPMPEVLYEEVLEVDERVVLHRGEAGTGTPVKGRTGDLLEVQQPVDLGALRGKLEGLL" \
                      "SRGIRSLAVVLMHSYTWAQHEQQVGVLARELGFTHVSLSSEAMPMVRIVPRGHTACADAY" \
                      "LTPAIQRYVQGFCRGFQGQLKDVQVLFMRSDGGLAPMDTFSGSSAVLSGPAGGVVGYSAT" \
                      "TYQQEGGQPVIGFDMGGTSTDVSRYAGEFEHVFEASTAGVTLQAPQLDINTVAAGGGSRL" \
                      "FFRSGLFVVGPESAGAHPGPACYRKGGPVTVTDANLVLGRLLPASFPCIFGPGENQPLSP" \
                      "EASRKALEAVATEVNSFLTNGPCPASPLSLEEVAMGFVRVANEAMCRPIRALTQARGHDP" \
                      "SAHVLACFGGAGGQHACAIARALGMDTVHIHRHSGLLSALGLALADVVHEAQEPCSLLYA" \
                      "PETFVQLDQRLSRLEEQCVDALQAQGFPRSQISTESFLHLRYQGTDCALMVSAHQHPATA" \
                      "RSPRAGDFGAAFVERYMREFGFVIPERPVVVDDVRVRGTGRSGLRLEDAPKAQTGPPRVD" \
                      "KMTQCYFEGGYQETPVYLLAELGYGHKLHGPCLIIDSNSTILVEPGCQAEVTKTGDICIS" \
                      "VGAEVPGTVGPQLDPIQLSIFSHRFMSIAEQMGRILQRTAISTNIKERLDFSCALFGPDG" \
                      "GLVSNAPHIPVHLGAMQETVQFQIQHLGADLHPGDVLLSNHPSAGGSHLPDLTVITPVFW" \
                      "PGQTRPVFYVASRGHHADIGGITPGSMPPHSTMLQQEGAVFLSFKLVQGGVFQEEAVTEA" \
                      "LRAPGKVPNCSGTRNLHDNLSDLRAQVAANQKGIQLVGELIGQYGLDVVQAYMGHIQANA" \
                      "ELAVRDMLRAFGTSRQARGLPLEVSSEDHMDDGSPIRLRVQISLSQGSAVFDFSGTGPEV" \
                      "FGNLNAPRAVTLSALIYCLRCLVGRDIPLNQGCLAPVRVVIPRGSILDPSPEAAVVGGNV" \
                      "LTSQRVVDVILGAFGACAASQGCMNNVTLGNAHMGYYETVAGGAGAGPSWHGRSGVHSHM" \
                      "TNTRITDPEILESRYPVILRRFELRRGSGGRGRFRGGDGVTRELLFREEALLSVLTERRA" \
                      "FRPYGLHGGEPGARGLNLLIRKNGRTVNLGGKTSVTVYPGDVFCLHTPGGGGYGDPEDPA" \
                      "PPPGSPPQALAFPEHGSVYEYRRAQEAV"

sequences['O14841']['ensembl_transcript_id'] = 'ENST00000618853'
sequences['O14841']['ensembl_protein_id']    = 'ENSP00000480476'

###################################################
def main():

	db, cursor = connect()
	switch_to_db(cursor, 'monogenic_development')
	for uniprot_id, updatable in sequences.iteritems():
		fixed_fields  = {'uniprot_id': uniprot_id}
		update_fields = updatable
		store_or_update (cursor, 'uniprot_seqs', fixed_fields, update_fields)
	return


#########################################
if __name__ == '__main__':
	main()
