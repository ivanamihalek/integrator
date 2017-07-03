#!/usr/bin/env python
# $Id: ncbiblast_urllib2.py 2673 2013-10-26 05:56:18Z hpm $
# ======================================================================
#
# Copyright 2009-2013 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ======================================================================
# NCBI BLAST (REST) Python client using urllib2 and
# xmltramp (http://www.aaronsw.com/2002/xmltramp/).
#
# Tested with:
#  Python 2.5.2 (Ubuntu 8.04 LTS)
#  Python 2.6.5 (Ubuntu 10.04 LTS)
#  Python 2.7.3 (Ubuntu 12.04 LTS)
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/sss/ncbi_blast_rest
# http://www.ebi.ac.uk/Tools/webservices/tutorials/python
# ======================================================================
# Base URL for service
baseUrl = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast'

# Load libraries
import platform, os, re, sys, time, urllib, urllib2, xmltramp2
from optparse import OptionParser

# Set interval for checking status
checkInterval = 10
# Output level
outputLevel = 1
# Debug level
debugLevel = 0
# Number of option arguments.
numOpts = len(sys.argv)

# Usage message
usage = "Usage: %prog [options...] [seqFile]"
description = """Rapid sequence database search programs utilizing the BLAST algorithm. For more information 
on NCBI BLAST refer to http://www.ebi.ac.uk/Tools/sss/ncbiblast"""
epilog = """For further information about the NCBI BLAST (SOAP) web service, see 
http://www.ebi.ac.uk/Tools/webservices/services/sss/ncbi_blast_soap."""
version = "$Id: ncbiblast_urllib2.py 2673 2013-10-26 05:56:18Z hpm $"
# Process command-line options
parser = OptionParser(usage=usage, description=description, epilog=epilog, version=version)
# Tool specific options
parser.add_option('-p', '--program', help='program to run')
parser.add_option('-D', '--database', help='database to search')
parser.add_option('--stype', default='protein', help='query sequence type')
parser.add_option('-m', '--matrix', help='scoring matrix')
parser.add_option('-E', '--exp', help='E-value threshold')
parser.add_option('-f', '--filter', action="store_true", help='low complexity sequence filter')
parser.add_option('-n', '--alignments', type='int', help='maximum number of alignments')
parser.add_option('-s', '--scores', type='int', help='maximum number of scores')
parser.add_option('-d', '--dropoff', type='int', help='dropoff score')
parser.add_option('--match_score', help='match/missmatch score')
parser.add_option('-o', '--gapopen', type='int', help='open gap penalty')
parser.add_option('-x', '--gapext', type='int', help='extend gap penalty')
parser.add_option('-g', '--gapalign', action="store_true", help='optimise gap alignments')
parser.add_option('--compstats', help='compositional adjustment/statistics mode')
parser.add_option('--seqrange', help='region within input to use as query')
parser.add_option('--sequence', help='input sequence file name')
# General options
parser.add_option('--email', help='e-mail address')
parser.add_option('--title', help='job title')
parser.add_option('--outfile', help='file name for results')
parser.add_option('--outformat', help='output format for results')
parser.add_option('--async', action='store_true', help='asynchronous mode')
parser.add_option('--jobid', help='job identifier')
parser.add_option('--polljob', action="store_true", help='get job result')
parser.add_option('--status', action="store_true", help='get job status')
parser.add_option('--resultTypes', action='store_true', help='get result types')
parser.add_option('--params', action='store_true', help='list input parameters')
parser.add_option('--paramDetail', help='get details for parameter')
parser.add_option('--quiet', action='store_true', help='decrease output level')
parser.add_option('--verbose', action='store_true', help='increase output level')
parser.add_option('--baseURL', default=baseUrl, help='Base URL for service')
parser.add_option('--debugLevel', type='int', default=debugLevel, help='debug output level')

(options, args) = parser.parse_args()

print options
print "*******************"
print args

# Increase output level
if options.verbose:
	outputLevel += 1

# Decrease output level
if options.quiet:
	outputLevel -= 1

# Debug level
if options.debugLevel:
	debugLevel = options.debugLevel


# Debug print
def printDebugMessage(functionName, message, level):
	if (level <= debugLevel):
		print >> sys.stderr, '[' + functionName + '] ' + message


# User-agent for request (see RFC2616).
def getUserAgent():
	printDebugMessage('getUserAgent', 'Begin', 11)
	# Agent string for urllib2 library.
	urllib_agent = 'Python-urllib/%s' % urllib2.__version__
	clientRevision = '$Revision: 2673 $'
	clientVersion = '0'
	if len(clientRevision) > 11:
		clientVersion = clientRevision[11:-2]
	# Prepend client specific agent string.
	user_agent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
		clientVersion, os.path.basename(__file__),
		platform.python_version(), platform.system(),
		urllib_agent
	)
	printDebugMessage('getUserAgent', 'user_agent: ' + user_agent, 12)
	printDebugMessage('getUserAgent', 'End', 11)
	return user_agent


# Wrapper for a REST (HTTP GET) request
def restRequest(url):
	printDebugMessage('restRequest', 'Begin', 11)
	printDebugMessage('restRequest', 'url: ' + url, 11)
	try:
		# Set the User-agent.
		user_agent = getUserAgent()
		http_headers = {'User-Agent': user_agent}
		req = urllib2.Request(url, None, http_headers)
		# Make the request (HTTP GET).
		reqH = urllib2.urlopen(req)
		result = reqH.read()
		reqH.close()
	# Errors are indicated by HTTP status codes.
	except urllib2.HTTPError, ex:
		# Trap exception and output the document to get error message.
		print >> sys.stderr, ex.read()
		raise
	printDebugMessage('restRequest', 'End', 11)
	return result


# Get input parameters list
def serviceGetParameters():
	printDebugMessage('serviceGetParameters', 'Begin', 1)
	requestUrl = baseUrl + '/parameters'
	printDebugMessage('serviceGetParameters', 'requestUrl: ' + requestUrl, 2)
	xmlDoc = restRequest(requestUrl)
	doc = xmltramp2.parse(xmlDoc)
	printDebugMessage('serviceGetParameters', 'End', 1)
	return doc['id':]


# Print list of parameters
def printGetParameters():
	printDebugMessage('printGetParameters', 'Begin', 1)
	idList = serviceGetParameters()
	for id in idList:
		print id
	printDebugMessage('printGetParameters', 'End', 1)


# Get input parameter information
def serviceGetParameterDetails(paramName):
	printDebugMessage('serviceGetParameterDetails', 'Begin', 1)
	printDebugMessage('serviceGetParameterDetails', 'paramName: ' + paramName, 2)
	requestUrl = baseUrl + '/parameterdetails/' + paramName
	printDebugMessage('serviceGetParameterDetails', 'requestUrl: ' + requestUrl, 2)
	xmlDoc = restRequest(requestUrl)
	doc = xmltramp2.parse(xmlDoc)
	printDebugMessage('serviceGetParameterDetails', 'End', 1)
	return doc


# Print description of a parameter
def printGetParameterDetails(paramName):
	printDebugMessage('printGetParameterDetails', 'Begin', 1)
	doc = serviceGetParameterDetails(paramName)
	print str(doc.name) + "\t" + str(doc.type)
	print doc.description
	for value in doc.values:
		print value.value,
		if str(value.defaultValue) == 'true':
			print 'default',
		print
		print "\t" + str(value.label)
		if (hasattr(value, 'properties')):
			for wsProperty in value.properties:
				print  "\t" + str(wsProperty.key) + "\t" + str(wsProperty.value)
	# print doc
	printDebugMessage('printGetParameterDetails', 'End', 1)


# Submit job
def serviceRun(email, title, params):
	printDebugMessage('serviceRun', 'Begin', 1)
	# Insert e-mail and title into params
	params['email'] = email
	if title:
		params['title'] = title
	requestUrl = baseUrl + '/run/'
	printDebugMessage('serviceRun', 'requestUrl: ' + requestUrl, 2)
	# Database requires special handling, so extract from params
	databaseList = params['database']
	del params['database']
	# Build the database data options
	databaseData = ''
	for db in databaseList:
		databaseData += '&database=' + db
	# Get the data for the other options
	requestData = urllib.urlencode(params)
	# Concatenate the two parts.
	requestData += databaseData
	printDebugMessage('serviceRun', 'requestData: ' + requestData, 2)
	# Errors are indicated by HTTP status codes.
	try:
		# Set the HTTP User-agent.
		user_agent = getUserAgent()
		http_headers = {'User-Agent': user_agent}
		req = urllib2.Request(requestUrl, None, http_headers)
		# Make the submission (HTTP POST).
		reqH = urllib2.urlopen(req, requestData)
		jobId = reqH.read()
		reqH.close()
	except urllib2.HTTPError, ex:
		# Trap exception and output the document to get error message.
		print >> sys.stderr, ex.read()
		raise
	printDebugMessage('serviceRun', 'jobId: ' + jobId, 2)
	printDebugMessage('serviceRun', 'End', 1)
	return jobId


# Get job status
def serviceGetStatus(jobId):
	printDebugMessage('serviceGetStatus', 'Begin', 1)
	printDebugMessage('serviceGetStatus', 'jobId: ' + jobId, 2)
	requestUrl = baseUrl + '/status/' + jobId
	printDebugMessage('serviceGetStatus', 'requestUrl: ' + requestUrl, 2)
	status = restRequest(requestUrl)
	printDebugMessage('serviceGetStatus', 'status: ' + status, 2)
	printDebugMessage('serviceGetStatus', 'End', 1)
	return status


# Print the status of a job
def printGetStatus(jobId):
	printDebugMessage('printGetStatus', 'Begin', 1)
	status = serviceGetStatus(jobId)
	print status
	printDebugMessage('printGetStatus', 'End', 1)


# Get available result types for job
def serviceGetResultTypes(jobId):
	printDebugMessage('serviceGetResultTypes', 'Begin', 1)
	printDebugMessage('serviceGetResultTypes', 'jobId: ' + jobId, 2)
	requestUrl = baseUrl + '/resulttypes/' + jobId
	printDebugMessage('serviceGetResultTypes', 'requestUrl: ' + requestUrl, 2)
	xmlDoc = restRequest(requestUrl)
	doc = xmltramp2.parse(xmlDoc)
	printDebugMessage('serviceGetResultTypes', 'End', 1)
	return doc['type':]


# Print list of available result types for a job.
def printGetResultTypes(jobId):
	printDebugMessage('printGetResultTypes', 'Begin', 1)
	resultTypeList = serviceGetResultTypes(jobId)
	for resultType in resultTypeList:
		print resultType['identifier']
		if (hasattr(resultType, 'label')):
			print "\t", resultType['label']
		if (hasattr(resultType, 'description')):
			print "\t", resultType['description']
		if (hasattr(resultType, 'mediaType')):
			print "\t", resultType['mediaType']
		if (hasattr(resultType, 'fileSuffix')):
			print "\t", resultType['fileSuffix']
	printDebugMessage('printGetResultTypes', 'End', 1)


# Get result
def serviceGetResult(jobId, type):
	printDebugMessage('serviceGetResult', 'Begin', 1)
	printDebugMessage('serviceGetResult', 'jobId: ' + jobId, 2)
	printDebugMessage('serviceGetResult', 'type: ' + type, 2)
	requestUrl = baseUrl + '/result/' + jobId + '/' + type
	result = restRequest(requestUrl)
	printDebugMessage('serviceGetResult', 'End', 1)
	return result


# Client-side poll
def clientPoll(jobId):
	printDebugMessage('clientPoll', 'Begin', 1)
	result = 'PENDING'
	while result == 'RUNNING' or result == 'PENDING':
		result = serviceGetStatus(jobId)
		print >> sys.stderr, result
		if result == 'RUNNING' or result == 'PENDING':
			time.sleep(checkInterval)
	printDebugMessage('clientPoll', 'End', 1)


# Get result for a jobid
def getResult(jobId):
	printDebugMessage('getResult', 'Begin', 1)
	printDebugMessage('getResult', 'jobId: ' + jobId, 1)
	# Check status and wait if necessary
	clientPoll(jobId)
	# Get available result types
	resultTypes = serviceGetResultTypes(jobId)
	for resultType in resultTypes:
		# Derive the filename for the result
		if options.outfile:
			filename = options.outfile + '.' + str(resultType['identifier']) + '.' + str(resultType['fileSuffix'])
		else:
			filename = jobId + '.' + str(resultType['identifier']) + '.' + str(resultType['fileSuffix'])
		# Write a result file
		if not options.outformat or options.outformat == str(resultType['identifier']):
			# Get the result
			result = serviceGetResult(jobId, str(resultType['identifier']))
			fh = open(filename, 'w');
			fh.write(result)
			fh.close()
			print filename
	printDebugMessage('getResult', 'End', 1)


# Read a file
def readFile(filename):
	printDebugMessage('readFile', 'Begin', 1)
	fh = open(filename, 'r')
	data = fh.read()
	fh.close()
	printDebugMessage('readFile', 'End', 1)
	return data

print "still alive"
# No options... print help.
if numOpts < 2:
	parser.print_help()
# List parameters
elif options.params:
	printGetParameters()
# Get parameter details
elif options.paramDetail:
	printGetParameterDetails(options.paramDetail)
# Submit job
elif options.email and not options.jobid:
	params = {}
	if len(args) > 0:
		if os.access(args[0], os.R_OK):  # Read file into content
			params['sequence'] = readFile(args[0])
		else:  # Argument is a sequence id
			params['sequence'] = args[0]
	elif options.sequence:  # Specified via option
		if os.access(options.sequence, os.R_OK):  # Read file into content
			params['sequence'] = readFile(options.sequence)
		else:  # Argument is a sequence id
			params['sequence'] = options.sequence
	# Booleans need to be represented as 1/0 rather than True/False
	if options.gapalign is not None:
		if options.gapalign:
			params['gapalign'] = True
		else:
			params['gapalign'] = False
	# Add the other options (if defined)
	if options.program:
		params['program'] = options.program
	if options.database:
		params['database'] = re.split('[ \t\n,;]+', options.database)
	if options.stype:
		params['stype'] = options.stype
	if options.matrix:
		params['matrix'] = options.matrix
	if options.exp:
		params['exp'] = options.exp
	if options.filter:
		params['filter'] = options.filter
	if options.alignments:
		params['alignments'] = options.alignments
	if options.scores:
		params['scores'] = options.scores
	if options.dropoff:
		params['dropoff'] = options.dropoff
	if options.match_score:
		params['match_score'] = options.match_score
	if options.gapopen:
		params['gapopen'] = options.gapopen
	if options.gapext:
		params['gapext'] = options.gapext
	if options.compstats:
		params['compstats'] = options.compstats

	# Submit the job
	jobid = serviceRun(options.email, options.title, params)
	if options.async:  # Async mode
		print jobid
	else:  # Sync mode
		print >> sys.stderr, jobid
		time.sleep(5)
		getResult(jobid)
# Get job status
elif options.status and options.jobid:
	printGetStatus(options.jobid)
# List result types for job
elif options.resultTypes and options.jobid:
	printGetResultTypes(options.jobid)
# Get results for job
elif options.polljob and options.jobid:
	getResult(options.jobid)
else:
	print >> sys.stderr, 'Error: unrecognised argument combination'
	parser.print_help()
