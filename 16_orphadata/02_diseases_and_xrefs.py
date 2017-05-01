#!/usr/bin/python

import xmltodict

def recursive_print(mistery_object):
	if not  mistery_object: return
	if type(mistery_object) is unicode or type(mistery_object) is str:
		print "is unicode or  string:"
		print mistery_object
	elif type(mistery_object) is list:
		print "is list:"
		for item in mistery_object:
			print   recursive_print(item)
	else:
		try:
			for k2, v2 in mistery_object.iteritems():
				if not v2: continue
				print
				print "key:", k2
				print recursive_print (v2)
		except:
			print "unrecognized type", mistery_object


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
		recursive_print (disorder)
		print
		exit()


#########################################
if __name__ == '__main__':
	main()
