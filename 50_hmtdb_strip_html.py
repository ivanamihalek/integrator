#!/usr/bin/python


from urllib2 import urlopen
from bs4 import BeautifulSoup


####################################
# this is not really working
def get_url_contents(url):
	ret = urlopen(url).read().decode('utf8')
	return  BeautifulSoup(ret,"lxml").get_text()


####################################
def main():
	# Write a utility function that takes a URL as its argument,
	# and returns the contents of the URL, with all HTML markup removed.
	# Use from urllib import request and then (this is antiquated)
	# request.urlopen('http://nltk.org/').read().decode('utf8') to access the contents of the URL.
	#html_text =  get_url_contents('http://www.mitomap.org/foswiki/bin/view/MITOMAP/PolymorphismsCoding')
	infile = open('/databases/mitomap/MITOMAP_PolymorphismsCoding.html')
	reading=False
	ct = 0
	for line in infile:
		print line
		if line [:len('"data":[')] == '"data":[':
			reading=True
		elif reading:
			print line
			ct += 1
			if ct==10: break
	infile.close()

	return True


####################################
if __name__ == '__main__':
	main()
