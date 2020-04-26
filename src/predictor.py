from tissuetypes import gettissues
from interactionconnector import expandnodes
from dataparser import Converter
import networkx as nx
import datetime

class Prediction:
	def __init__(self, query):
		self.query = query
		self.network = {}
		self.predictions = {}
		self.analysis = {}
		#self.network contains tissue to network pairs
		#self.predictions contains tissue to list of list of nodes (list of nodes = path) pairs
		#self.analysis contains stats to value pairs
		self.status = 0
		#self.status values:
			#0 = Initialized, pending
			#1 = Complete prediction with no errors
			#2 = Complete analysis with no errors
			#-2 = Some protein not found in database
			#-3 = No compatible common tissue types
			#-1 = Unspecified error(reserved value)
		
	def predict(self, _logger, _converter, _stringdb, _padb):
	
		_logger.log("Predicting pathway between {0} and {1}".format(self.query[0], self.query[1]), toprint = True, timestamp = True)
	
		#parse query
		startacc, endacc, iterations = self.query

		#convert to ACC if necessary from ENSG/ENSP
		if startacc[0] == "E":
			inittype = ""
			if startacc[3] == "G":
				inittype = "ENSEMBL_ID"
				
			if startacc[3] == "P":
				inittype = "ENSEMBL_PRO_ID"
				
			d = _converter.mass_convert(_logger, [startacc, endacc], inittype, "ACC")
			startacc = d[startacc]
			endacc = d[endacc]


		#initialize network
		_initnodes = set([startacc, endacc])
		_initedges = set()
		starttissues = gettissues(_logger, startacc, _padb)
		endtissues = gettissues(_logger, endacc, _padb)
		if starttissues == set() or endtissues == set():
			_logger.log("Some proteins not found in tissuetypes database\n{0}: {1}\n{2}: {3}".format(startacc, starttissues, endacc, endtissues), toprint = True)
			self.network = {}
			self.predictions = {}
			self.status = -2
			
			return
			
		_inittissues = starttissues & endtissues
		
		#verify viability
		if len(_inittissues) < 1:
			_logger.log("\nNo common tissues between {0} and {1}\n{0}: {2}\n{1}: {3}\n"
				.format(startacc, endacc, starttissues, endtissues), toprint = True)	

			self.network = {}
			self.predictions = {}
			self.status = -3
				
			return
		
		#iterate network
		for t in _inittissues: 
		
			self.predictions[t] = []
		
			#init network
			_nodes = set(_initnodes)
			_edges = set(_initedges)
			
			#iterate
			for i in range(iterations):
				expandnodes(_logger, _converter, _nodes, _edges, _stringdb, _padb, t, i)
				
			#graph creation
			G = nx.Graph()
			G.add_nodes_from(_nodes)
			G.add_edges_from(_edges)
			
			self.network[t] = G
			
			#solve
			try:
				self.predictions[t] = list(nx.all_shortest_paths(G, source = startacc, target = endacc))
			except nx.exception.NetworkXNoPath as e:
				self.predictions[t] = []
			
		self.status = 1
			
	def analyze(self, _logger, analyzer, *args):
	
		if self.status != 1:
		
			_logger.log("Error! Status = {0}".format(self.status), toprint = True, timestamp = True)
		
			return
	
		_logger.log("Analyzing predicted pathways between {0} and {1}".format(self.query[0], self.query[1]), toprint = True, timestamp = True)
		
		stats = analyzer(_logger, self, *args)
		self.analysis = stats
		self.status = 2
	
		return stats
		
	def tostring(self, _genenames):
	
		if self.status != 2:
			return "Error! Status = {0}".format(self.status)
	
		strpredictions = "--- Predicted pathways between {0} and {1} ---\n\n".format(self.query[0], self.query[1])
		
		for t in self.predictions:
		
			tempstr = "Pathways in " + t + ":"
			count = 0
			
			for path in self.predictions[t]:
				
				count += 1
				tempstr += "\n >> #{0}:\t{1}\t".format(count, "\t".join(path))
				
				for protein in path:
					tempstr += "\t{0}".format(_genenames[protein])
				
			strpredictions += tempstr + "\n"
			
		stranalysis = "--- Prediction Analysis Results ---\n\n"
		
		for stat in self.analysis:
			val = self.analysis[stat]
		
			if isinstance(val, dict):
			
				tempstr = " > {0}:".format(stat)
				
				for protein in val:
					data = val[protein]
					if data == set():
						data = ""
					
					tempstr += "\n    @ {0}: {1}".format(protein, data)
		
			else:
				if val == set():
					val = ""
				tempstr = " > {0}: {1}".format(stat, val)
				
			stranalysis += tempstr + "\n"
			
		return strpredictions + "\n\n" + stranalysis


def predict(_logger, _converter, _stringdb, _padb, _genenames, query, analyzer, *args):

	_logger.log("Generating Predictor Object for query {0}".format(query), toprint = True, timestamp = True)
	
	P = Prediction(query)
	P.predict(_logger, _converter, _stringdb, _padb)
	P.analyze(_logger, analyzer, *args)
	
	return P
	