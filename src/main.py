import sys
import datetime
import csv
from predictor import predict
from dataparser import databasesinit, save, Converter
from analyzer import analyze, mass_analyze
from Logger import Logger

#return help message on 0, args required on 1
#edit reqargs as programme is built
def help(i):

	reqargs = [
		"main.py", 
		"queryfile", #queryfile contains expectedpathways
		"n_iter", 
		"stringdb_threshold", 
		"proteinatlasdb_threshold",
		]
		
	str = " ".join(reqargs)
	result = (str, len(reqargs))
	return result[i]

def main(args):

	#process args
	queryfile, iterations, s_threshold, p_threshold = args
	iterations = int(iterations)
	s_threshold = float(s_threshold)
	p_threshold = float(p_threshold)

	#logger and converter initialization
	_logger = Logger(str(datetime.datetime.now()))
	_converter = Converter()

	#parse queryfile into list of queries
	queries = []
	with open(queryfile, "r") as f:
		reader = csv.reader(f)
		for row in reader:
			queries.append(row)
			
	#convert all to acc		
	tocheck = []
	for expectedpath in queries:
		tocheck += expectedpath
	
	#check code type of input
	code = queries[0][0]
	if code[0] == "E":
		inittype = ""
		if code[3] == "G":
			inittype = "ENSEMBL_ID"
			
		if code[3] == "P":
			inittype = "ENSEMBL_PRO_ID"
	
	#run all together to minimize API pings
	converted = _converter.mass_convert(_logger, tocheck, inittype, "ACC")
	
	#then convert
	for i in range(len(queries)):
		for j in range(len(queries[i])):
			queries[i][j] = converted[queries[i][j]]

	#databases initialization
	_stringdb, _proteinatlas, _genenames = databasesinit(_logger, s_threshold, p_threshold)
	
	#output initialization
	_output = {}

	#run prediction
	for q in queries:
		query = (q[0], q[-1], iterations)
		_output[query] = predict(_logger, _converter, _stringdb, _proteinatlas, _genenames, query, analyze, q)
	
	#analyze en-masse
	collated, printablecollated, totalanalysis = mass_analyze(_logger, _output.values(), _genenames)
	
	#output analysis data, then save _logger output
	save(_logger, _output, printablecollated, totalanalysis, s_threshold, p_threshold, _genenames)


if __name__ == "__main__":
	if len(sys.argv) != help(1):
		print(help(0))
	else:
		main(sys.argv[1:])