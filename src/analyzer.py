import networkx as nx
import itertools

def analyze(_logger, Predictor, *args):

	expectedpath = args[0]

	stats = {
		"length": set(),
		"n_tissues": -1,
		"largest tissue pathway set": -1,
		"smallest tissue pathway set": -1,
		"FP <= E": -1,
		"expected path present": [],
		"expected protein presence": {},
		"expected protein position in pathway": {}
		}
		
	
		
	#compute length
	for t in Predictor.predictions:
		if len(Predictor.predictions[t]) > 0:
			stats["length"].add(len(Predictor.predictions[t][0])-1)
	stats["length"] = list(stats["length"])
	
	#compute n_tissues
	stats["n_tissues"] = len(Predictor.predictions)
	
	#compute largest and smallest tissue pathway set
	templs = []
	for t in Predictor.predictions:
		templs.append(len(Predictor.predictions[t]))
	stats["largest tissue pathway set"] = max(templs)
	stats["smallest tissue pathway set"] = min(templs)
	
	#compute false positives of length <= expected path
	tempsum = 0
	for t in Predictor.network:
		tempsum += len(list(nx.all_simple_paths(Predictor.network[t], source = Predictor.query[0], target = Predictor.query[1], cutoff = 3))) 
	stats["FP <= E"] = tempsum
	
	#compute expected path presence
	expectededgelist = []
	for i in range(len(expectedpath)-1):
		expectededgelist.append((expectedpath[i], expectedpath[i+1]))
		
	presentin = []
	for t in Predictor.network:
		G = Predictor.network[t]
		networkedgelist = list(G.edges)
		
		presence = True
		for edge in expectededgelist:
			if edge not in networkedgelist:
				presence = False
				break
		
		if presence:
			presentin.append(t)
			
	if presentin == []:
		presentin = False
			
	stats["expected path present"] = presentin
	
	#compute expected protein presence and protein position in pathway
	for protein in expectedpath:
		tissuespresent = set()
		positions = set()
		
		for t in Predictor.predictions:
			for path in Predictor.predictions[t]:
				if protein in path:
					tissuespresent.add(t)
					positions.add(path.index(protein))
					
		stats["expected protein presence"][protein] = tissuespresent
		stats["expected protein position in pathway"][protein] = positions
		
	#compute expected protein presence in combinations:
	for i in range(2, len(expectedpath) + 1):
		combis = itertools.combinations(expectedpath, i)
		for combi in combis:
			intersection = stats["expected protein presence"][combi[0]]
			for protein in combi[1:]:
				intersection = intersection & stats["expected protein presence"][protein]
				
			stats["expected protein presence"][combi] = intersection
		
	return stats
	
def mass_analyze(_logger, predictions, _genenames):
	
	_logger.log("Conducting Meta-Analysis", toprint=True, timestamp=True)
	
	#the first part aggregates all the stats
	rows = []
	printout = []
	headers = ["Start", "End", "Status"]
	for prediction in predictions:
	
		if prediction.status != 2:
			currrow = [prediction.query[0], prediction.query[1], prediction.status]
			rows.append(currrow)
			for i in range(len(currrow)):
				currrow[i] = str(currrow[i])
			printout.append("\t".join(currrow + [prediction.tostring(_genenames)]))
		
		else:
			if len(headers) == 3:
				for key in prediction.analysis:
					if not isinstance(prediction.analysis[key], dict):
						headers.append(key)
				
				rows = headers + rows
				printout = ["\t".join(headers)] + printout
				
			currrow = [prediction.query[0], prediction.query[1], prediction.status]
			for stat in headers[3:]:
				currrow.append(prediction.analysis[stat])
				
			rows.append(currrow)
			for i in range(len(currrow)):
				currrow[i] = str(currrow[i])
			printout.append("\t".join(currrow))
	
	#the second part does some simple meta analysis
	metastats = {
		"average_length": 0,
		}
		
	#calculate average_len:
	n = 0
	validpredictions = list(filter(lambda p: p.status == 2, predictions))
	for prediction in validpredictions:
		lengths = prediction.analysis["length"]
		if len(lengths) > 0:
			n += sum(lengths)/len(lengths)
		
	if len(validpredictions) > 0:
		metastats["average_length"] = n/len(validpredictions)
		
	return rows, printout, metastats
	