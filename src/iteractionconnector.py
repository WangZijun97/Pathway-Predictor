from dataparser import mass_convert
from tissuetypes import gettissues
import datetime

def expandnodes(_logger, _nodes, _edges, _stringdb, _padb, _tissue, iter):
	_logger.log("[{0}] Expanding Nodes | Iteration {1}".format(datetime.datetime.now(), iter), toprint = True)
	
	#convert all nodes to STRING codes, stored in tempnodes
	lsnodes = list(_nodes)
	tempnodes = set()
	tempnodes.update(list(mass_convert(_logger, lsnodes, "ACC", "STRING_ID").values()))
	
	newnodes = set(tempnodes)
	newedges = set()
	
	#retrieving interactors and adding to newedges from stringdb
	for node in tempnodes:
		newnodes.update(_stringdb[node])
		for interactor in _stringdb[node]:
			newedges.add((node, interactor))
			
	
	
	#remove redundant edges
	newedges = list(set(newedges))
	
	#convert all nodes back to UniProt ACC
	convertion_table = mass_convert(_logger, list(newnodes), "STRING_ID", "ACC")
	_logger.log("[{0}] Conversion Query Complete. To check {1} nodes and {2} edges.".format(datetime.datetime.now(), len(newnodes), len(newedges)), toprint=True)
	convertededges = set()
	
	count = 0
	for edge in newedges:
		count += 1
	
		#filter if not found in convertion_table
		if edge[1] not in convertion_table:
			continue
	
		convertededges.add((convertion_table[edge[0]], convertion_table[edge[1]]))
		
		#debugging
		"""
		if count % 10000 == 0:
			_logger.log("[{0}] Converted {1} edges of {2}".format(datetime.datetime.now(), count, len(newedges)), toprint = True)
			"""
			
	convertednodes = set(list(convertion_table.values()))
	
	#filter based on tissuetype
	unwantededges = set()
	unwantednodes = set()
	
	_logger.log("[{0}] Conversion Complete. Filtering {1} edges".format(datetime.datetime.now(), len(convertededges)), toprint = True)
	
	count = 0
	for edge in convertededges:
		count += 1
		if _tissue not in gettissues(_logger, edge[1], _padb):
			unwantededges.add(edge)
			unwantednodes.add(edge[1])
		
		#debugging
		"""
		if count % 10000 == 0:
			_logger.log("[{0}] Filtered {1} edges of {2}".format(datetime.datetime.now(), count, len(convertededges)), toprint = True)
			"""
			
	convertednodes = convertednodes - unwantednodes
	convertededges = convertededges - unwantededges
	
	#update original containers
	_nodes.update(convertednodes)
	_edges.update(convertededges)
	
	_logger.log("[{0}] Iteration {1} Complete".format(datetime.datetime.now(), iter), toprint = True)	
	
	return _nodes, _edges
	