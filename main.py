import collections
import json
import networkx as nx
import requests
import sys
import time
from xmlparse import getName

#Open Memo tables to save time on repeat API queries
#Memoization table for retrieve()
with open("interactions.json", "r+") as fRetrieve:
	RETRIEVE_MEMOIZATION_TABLE = json.load(fRetrieve)

#Memoization table for getName()
with open("accname.json", "r+") as fAccName:
	ACCTONAME_MEMOIZATION_TABLE = json.load(fAccName)
	
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
	
#node expansion function
#nodelist is a set
#edgelist is a list
def expandNodes(nodelist, edgelist):
	
	newnodes = set()
	for node in nodelist:
		for entity in retrieve(node)["entities"]:
			for interactor in entity["interactors"]:
				newnodes.update([interactor["acc"]])
				edgelist += [(entity["acc"], interactor["acc"])]
				
	nodelist.update(newnodes)
	return nodelist, edgelist
	
def display(ls, front):
	str = "{0}\n--------------------\n\n".format(front)
	for acc in ls:
		str += " > {0} ({1})\n".format(getName(acc, ACCTONAME_MEMOIZATION_TABLE), acc)
	return str
	
def main(args):

	#process args
	startacc, endacc, iterations = args
	iterations = int(iterations)

	#node and edge storage
	_nodes = set()
	_edges = []
	
	#initialize
	_nodes.update([startacc, endacc])
	
	#iterate
	for i in range(iterations):
		expandNodes(_nodes, _edges)
		
	#graph creation
	G = nx.Graph()
	G.add_nodes_from(_nodes)
	G.add_edges_from(_edges)
	
	#solve
	try:
		count = 0
		for path in nx.all_shortest_paths(G, source = startacc, target = endacc):
			count += 1
			print(display(path, "Pathway #" + str(count)))
	except nx.exception.NetworkXNoPath as e:
		print(e)
	
	#print("Nodes: " + str(_nodes))
	#print("Edges: " + str(_edges))
	
	
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("main.py startacc endacc n_iterations")
	else:	
		main(sys.argv[1:])
		
	#close files
	with open("interactions.json", "w") as f:
		json.dump(RETRIEVE_MEMOIZATION_TABLE, f)
		f.close()
		
	with open("accname.json", "w") as f:
		json.dump(ACCTONAME_MEMOIZATION_TABLE, f)
		f.close()
	
	
	
	
	
	