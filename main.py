import collections
import json
import networkx as nx
import requests
import sys
import time
from datetime import datetime
from xmlparse import getName, getNameFromSTRING
from tissuecodes import gettissues

#Open Memo tables to save time on repeat API queries
#Memoization table for retrieve()
with open("interactions.json", "r+") as fRetrieve:
	RETRIEVE_MEMOIZATION_TABLE = json.load(fRetrieve)

#Memoization table for getName()
with open("accname.json", "r+") as fAccName:
	ACCTONAME_MEMOIZATION_TABLE = json.load(fAccName)
	
#Memoization table for getNameFromSTRING()
with open("stringid.json", "r+") as fStringId:
	STRINGTOACC_MEMOIZATION_TABLE = json.load(fStringId)
	
#Memoization table for retrieve using STRING db implementation
with open("interactions_string.json", "r+") as fRetrieveString:
	RETRIEVE_MEMOIZATION_TABLE_STRING = json.load(fRetrieveString)

#Memoization table for retrieving data from proteinatlas.org	
with open("proteinatlas.json", "r+") as fProteinAtlas:
	PROTEIN_ATLAS_TABLE = json.load(fProteinAtlas)
	
"""
#Reactome implementation
def retrieve(acc):
	url = detailedInfoUrl(acc)
	if url in RETRIEVE_MEMOIZATION_TABLE:
		return RETRIEVE_MEMOIZATION_TABLE[url]
	else:
		data = requests.get(url).json()
		RETRIEVE_MEMOIZATION_TABLE[url] = data
		return data
	
def detailedInfoUrl(acc):
	return "https://reactome.org/ContentService/interactors/static/molecule/{0}/details".format(acc)
"""

class Counter:
	def __init__(self, init):
		self.i = init
		
	def plus(self):
		self.i += 1
		
class Logger:
	def __init__(self, time):
		self.loglist = []
		self.name = str(time).replace(":",".")
		
	def log(self, str):
		self.loglist.append(str)
		return str
		
	def printout(self):
		return "\n".join(self.loglist)

#STRING-db implementation
def retrieve(_logger, id, counter):s
	
	url = urlBuilder(id)
	if url in RETRIEVE_MEMOIZATION_TABLE_STRING:
		return RETRIEVE_MEMOIZATION_TABLE_STRING[url]
	else:
		print(_logger.log("[{0}] Retrieving String Interactors of {1}".format(datetime.now(), id)))
		data = requests.get(url).json()
		if counter:
			counter.plus()
			print(_logger.log("HTML Retrieval #{0}: {1}".format(counter.i, url)))
		RETRIEVE_MEMOIZATION_TABLE_STRING[url] = data
		return data
	
def urlBuilder(id):
	return "https://string-db.org/api/json/interaction_partners?identifiers={0}".format(id)
	
#node expansion function
#nodelist is a set
#edgelist is a list
def expandNodes(_logger, nodelist, edgelist, counter, tissue, iter):
	print(_logger.log("[{0}] Expanding Nodes | Iteration {1}".format(datetime.now(), iter)))
	
	newnodes = set()
	"""
	#Reactome Implementation
	for node in nodelist:
		for entity in retrieve(node)["entities"]:
			for interactor in entity["interactors"]:
				newnodes.update([interactor["acc"]])
				edgelist += [(entity["acc"], interactor["acc"])]
				
	"""
	
	#STRING-db implementation
	for node in nodelist:
		for interactor in retrieve(_logger, node, counter):
			newnodes.update([interactor["stringId_B"]])
			#newnodes.update([interactor["stringId_A"]])
			edgelist += [(interactor["stringId_A"], interactor["stringId_B"])]
	
	nodelist.update(filterNodes(_logger, newnodes, tissue))
	return nodelist, edgelist
	
#newnode processing - check tissue area
def filterNodes(_logger, nodelist, tissuetype):
	print(_logger.log("[{0}] Filtering nodes | Tissue type: {1} | Nodes to check: {2}".format(datetime.now(), tissuetype, len(nodelist))))
	filtered = set()
	for i in nodelist:
		if tissuetype in gettissues(_logger, i, PROTEIN_ATLAS_TABLE):
			filtered.update([i])
	return filtered
	
def display(ls, front):
	str = "{0}\n--------------------\n\n".format(front)
	for id in ls: #change to acc for reactome
		"""
		#Reactome implementation
		str += " > {0} ({1})\n".format(getName(acc, ACCTONAME_MEMOIZATION_TABLE), acc)
		"""
		
		#String-db implementation
		#str += " > {0} ({1})\n".format(getNameFromSTRING(id, STRINGTOACC_MEMOIZATION_TABLE, ACCTONAME_MEMOIZATION_TABLE), id)
		#str += " > {0}\n".format(id)
		str += " > {0} ({1})\n".format(getNameFromSTRING(id, STRINGTOACC_MEMOIZATION_TABLE), id)
	return str
	
def main(args):

	#process args
	startacc, endacc, iterations, _logger = args
	iterations = int(iterations)

	#node and edge storage
	_nodes = set()
	_edges = []
	
	#initialize
	_nodes.update([startacc, endacc])
	_counter = Counter(0)
	_tissueset = gettissues(_logger, startacc, PROTEIN_ATLAS_TABLE) & gettissues(_logger, endacc, PROTEIN_ATLAS_TABLE)
	
	
	#check if viable (same tissues existant)
	if len(_tissueset) < 1:
		print(_logger.log("\nNo common tissues between {0} and {1}\n{0}: {2}\n{1}: {3}\n"
			.format(startacc, endacc, gettissues(_logger, startacc, PROTEIN_ATLAS_TABLE), gettissues(_logger, endacc, PROTEIN_ATLAS_TABLE))))
		return
	
	#conduct solver for each tissue type
	for _tissue in _tissueset:
	
		_nodesspecific = set(_nodes)
		_edgesspecific = list(_edges)
	
		#iterate
		for i in range(iterations):
			expandNodes(_logger, _nodesspecific,_edgesspecific, _counter, _tissue, i)
			
		#graph creation
		G = nx.Graph()
		G.add_nodes_from(_nodesspecific)
		G.add_edges_from(_edgesspecific)
		
		#solve
		print(_logger.log("\n---------- {0} ----------\n".format(_tissue)))
		try:
			count = 0
			for path in nx.all_shortest_paths(G, source = startacc, target = endacc):
				count += 1
				print(_logger.log(display(path, "Pathway #" + str(count))))
		except nx.exception.NetworkXNoPath as e:
			print(_logger.log(str(e)))
		
		#print("Nodes: " + str(_nodesspecific))
		#print("Edges: " + str(_edgesspecific))
	
	
if __name__ == "__main__":
	_logger = Logger(datetime.now())
	if len(sys.argv) != 4:
		print(_logger.log("main.py startacc endacc n_iterations"))
	else:
		print(_logger.log("[{0}] Files opened, running algo\n".format(datetime.now())))
		main(sys.argv[1:]+[_logger])
		
	#close files
	print(_logger.log("\n--------------------\n\n"))
	print(_logger.log("[{0}] Closing files".format(datetime.now())))
	
	with open("interactions.json", "w") as f:
		json.dump(RETRIEVE_MEMOIZATION_TABLE, f)
		f.close()
		
	with open("accname.json", "w") as f:
		json.dump(ACCTONAME_MEMOIZATION_TABLE, f)
		f.close()
		
	with open("stringid.json", "w") as f:
		json.dump(STRINGTOACC_MEMOIZATION_TABLE, f)
		f.close()
	
	with open("interactions_string.json", "w") as f:
		json.dump(RETRIEVE_MEMOIZATION_TABLE_STRING, f)
		f.close()
		
	with open("proteinatlas.json", "w") as f:
		json.dump(PROTEIN_ATLAS_TABLE, f)
		f.close()
		
		
	print(_logger.log("[{0}] Files Closed".format(datetime.now())))
	with open("./output/{0}.txt".format(_logger.name), "w") as f:
		f.write(_logger.printout())
		f.close()
	
	
	
	
	
	