#!/usr/bin/python

import xmltodict

def recursive_print(mistery_object, depth):
	if not  mistery_object: return
	if type(mistery_object) is unicode or type(mistery_object) is str:
		#print "\t"*depth, "is unicode or string:"
		print "\t"*depth, mistery_object
	elif type(mistery_object) is list:
		#print "\t"*depth, "is list:"
		for item in mistery_object:
			recursive_print(item, depth+1)
	else:
		try:
			for k2, v2 in mistery_object.iteritems():
				if not v2: continue
				print "\t"*depth, "key:", k2
				recursive_print (v2, depth + 1)
		except:
			print "\t"*depth, "unrecognized type", mistery_object


##########################################
def main():

	with open('/databases/orphadata/rare_diseases_and_xrefs.orpahdata.xml') as fd:
		doc = xmltodict.parse(fd.read())
	print doc.keys()
	print doc['JDBOR'].keys()
	print doc['JDBOR']['DisorderList'].keys()

	# doc['JDBOR']['DisorderList']['Disorder'] is a list
	disorders = doc['JDBOR']['DisorderList']['Disorder']
	for disorder in disorders:
		recursive_print (disorder, 0)
		print
		exit()


#########################################
if __name__ == '__main__':
	main()
