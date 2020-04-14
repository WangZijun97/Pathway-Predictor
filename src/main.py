from Logger import Logger
from dataparser import string_db_dict_gen, proteinatlas_db_dict_gen_PROTEIN, mass_convert, proteinatlas_db_dict_gen
from tissuetypes import gettissues
from iteractionconnector import expandnodes
import datetime
import networkx as nx
import sys

def display(ls, front, db):
	str = "{0}\n--------------------\n\n".format(front)
	for id in ls: 
		str += " > {0} ({1})\n".format(db[id], id) #other format is db[id]["Gene"]
	return str

def main(args):

	#process args
	startacc, endacc, iterations = args
	iterations = int(iterations)
	
	#init logger
	_logger = Logger(str(datetime.datetime.now()) + " ({0} to {1})".format(startacc, endacc))
	_logger.log("[{0}] Beginning Query: Predicting Pathway from {1} to {2} in {3} iters".format(datetime.datetime.now(), startacc, endacc, iterations), toprint = True)
	
	#initial convertion
	if startacc[0] == "E":
		inittype = ""
		if startacc[3] == "G":
			inittype = "ENSEMBL_ID"
			
		if startacc[3] == "P":
			inittype = "ENSEMBL_PRO_ID"
			
		d = mass_convert(_logger, [startacc, endacc], inittype, "ACC")
		startacc = d[startacc]
		endacc = d[endacc]
		
	
	#initializations
	_initnodes = set([startacc, endacc])
	_initedges = set()
	_stringdb = string_db_dict_gen(_logger, 900)
	_proteinatlas, _genenames = proteinatlas_db_dict_gen(_logger)
	_inittissues = gettissues(_logger, startacc, _proteinatlas) & gettissues(_logger, endacc, _proteinatlas)
	_pathways = {}
	
	_logger.log("[{0}] Initial tissues: {1}".format(datetime.datetime.now(), _inittissues), toprint = True)
		
	#verify if viable (same tissues existant)
	if len(_inittissues) < 1:
		_logger.log("\nNo common tissues between {0} and {1}\n{0}: {2}\n{1}: {3}\n"
			.format(startacc, endacc, gettissues(_logger, startacc, _proteinatlas), gettissues(_logger, endacc, _proteinatlas)), toprint = True)
			
		#close stuff up
		_logger.log("[{0}] End of Query".format(datetime.datetime.now()), toprint = True)
		with open("../output/{0}.txt".format(_logger.name), "w") as f:
			f.write(_logger.printout())
			f.close()
		return
		
	#conduct network solving for each tissue type
	for _tissue in _inittissues: 
	
		_pathways[_tissue] = []
	
		#init network
		_nodes = set(_initnodes)
		_edges = set(_initedges)
		
		#iterate
		for i in range(iterations):
			expandnodes(_logger, _nodes, _edges, _stringdb, _proteinatlas, _tissue, i)
			
		#graph creation
		G = nx.Graph()
		G.add_nodes_from(_nodes)
		G.add_edges_from(_edges)
		
		#solve
		_logger.log("\n------ {0} ------\n".format(_tissue), toprint = True)
		try:
			count = 0
			for path in nx.all_shortest_paths(G, source = startacc, target = endacc):
				count += 1
				result = display(path, "Pathway #" + str(count), _genenames)
				_pathways[_tissue].append(result)
				_logger.log(result, toprint = True)
		except nx.exception.NetworkXNoPath as e:
			_logger.log(str(e) + "\n", toprint = True)
			
	#reprint all results
	for _tissue in _pathways:
		for pathway in _pathways[_tissue]:
			_logger.log(">>> {0} <<<\n\n{1}".format(_tissue, pathway), toprint = True)
			
	#save logger
	_logger.log("[{0}] End of Query".format(datetime.datetime.now()), toprint = True)
	with open("../output/{0}.txt".format(_logger.name), "w") as f:
		f.write(_logger.printout())
		f.close()
		
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("main.py startacc endacc n_iter")
	else:
		main(sys.argv[1:])
		
		
		
	