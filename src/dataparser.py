import csv
import json
import requests
import datetime
from Logger import Logger #only for testing

def string_db_dict_gen(_logger, confidencelimit):
	_logger.log("[{0}] Opening String-DB".format(datetime.datetime.now()), toprint = True)
	
	d = {}
	
	with open("../data/9606.protein.links.v11.0.txt", "r+") as f:
		reader = csv.reader(f, delimiter = " ")
		for row in reader:
			if row[2] == "combined_score":
				continue
			if float(row[2]) < confidencelimit:
				continue
			if row[0] not in d:
				d[row[0]] = set()
			if row[1] not in d:
				d[row[1]] = set()
			d[row[0]].add(row[1])
			d[row[1]].add(row[0])
			
	return d
	
def test_string_db_dict_gen(): #passed
	start = datetime.datetime.now()
	d = string_db_dict_gen()
	for i in sorted(d.keys()):
		print(i, d[i])
		print("\n\n")
	print(datetime.datetime.now() - start)
		
def proteinatlas_db_dict_gen(_logger):
	_logger.log("[{0}] Opening ProteinAtlas-DB".format(datetime.datetime.now()), toprint = True)
	
	d = {}
	names = {}
	
	with open("../data/proteinatlas.json", "r+") as f:
		jsondata = json.load(f)
	
	for i in jsondata:
		for code in i["Uniprot"]:
			d[code] = i
			names[code] = i["Gene"]
			
	return d, names

#this version uses pa-db protein expression data, and outputs already processed code to tissue list pairs
def proteinatlas_db_dict_gen_PROTEIN(_logger):
	_logger.log("[{0}] Opening ProteinAtlas-DB (Protein Expression)".format(datetime.datetime.now()), toprint = True)
	
	d = {}
	n = {}
	
	with open("../data/normal_tissue.tsv", "r") as f:
		reader = csv.reader(f, delimiter = "\t")
		for row in reader:
			if row[0] == "Gene":
				continue
			if row[0] not in d:
				d[row[0]] = []
			if row[0] not in n:
				n[row[0]] = row[1]
			if (row[4] == "High") and (row[5] != "Uncertain"):
				d[row[0]].append(row[2])
				
	keys = list(d.keys())
	result = {}
	name_convert = {}
				
	#convert keys to uniprot
	conversion_table = mass_convert(_logger, keys, "ENSEMBL_ID", "ACC")
	for i in d:
		if i in conversion_table:
			result[conversion_table[i]] = d[i]
			name_convert[conversion_table[i]] = n[i]
				
	return result, name_convert
		
def test_proteinatlas_db_dict_gen(): #passed, 18875 entries with uniprot codes
	d = proteinatlas_db_dict_gen(1)
	for i in sorted(d.keys()):
		print(d[i]["Gene"])
	print(len(d.keys()))
	
def test_proteinatlas_db_dict_gen_PROTEIN(): #passed, only 15313 entries available
	d = proteinatlas_db_dict_gen_PROTEIN(Logger(1))
	
	for i in d:
		print(i, d[i])
	print(len(d.keys()))
	


def mass_convert(_logger, ls_codes, type1, type2):
	_logger.log("[{0}] Converting {1} Codes of {2} to {3}".format(datetime.datetime.now(), len(ls_codes), type1, type2), toprint = True)
	
	#open memoi
	with open("../data/conversions.json", "r+") as f:
		conversions = json.load(f)
		f.close()
		
	d = {}
	toquery = []
	
	#check if entry is present for conversion type in memoi
	#initialize if not, and skip extraction
	if type1 not in conversions:
		conversions[type1] = {}
		toquery = ls_codes
		
	if type2 not in conversions[type1]:
		conversions[type1][type2] = {}
		toquery = ls_codes
		
	#extract if present
	else:
			
		for code in ls_codes:
			if code in conversions[type1][type2]:
				d[code] = conversions[type1][type2][code]
			else:
				toquery.append(code)
				
	#run uniprot query
	if toquery != []:
		queried = mass_convert_helper(_logger, toquery, type1, type2, recursed = False)
		d.update(queried)
		conversions[type1][type2].update(queried)
	
	
		#update memoi file
		with open("../data/conversions.json", "w") as f:
			json.dump(conversions, f)
			f.close()
	
	return d
	
	
def mass_convert_helper(_logger, ls_codes, type1, type2, recursed = False):
	_logger.log("[{0}] Converting {1} Codes of {2} to {3} via Uniprot".format(datetime.datetime.now(), len(ls_codes), type1, type2), toprint = not recursed)
	_logger.log(" > " + str(ls_codes))
	
	d = {}
	
	if len(ls_codes) > 200:
		for i in range(len(ls_codes)//200):
			d.update(mass_convert_helper(_logger, ls_codes[i*200:(i+1)*200], type1, type2, recursed = True))
		ls_codes = ls_codes[200 * (len(ls_codes)//200):]
		
	
	urlhead = "https://www.uniprot.org/uploadlists/?"
	urlparams = "from={0}&to={1}&format=tab&query={2}".format(type1, type2, " ".join(ls_codes))
	url = urlhead + urlparams
	_logger.log("[{0}] Querying uniprot.org with {1} codes".format(datetime.datetime.now(), len(ls_codes)), toprint = True)
	text = requests.get(url).text
	#print(text)
	ls = text.split("\n")
	for i in range(1, len(ls)):
		if not ls[i]:
			continue
		ls[i] = ls[i].split("\t")
		d[ls[i][0]] = ls[i][1]
		
	return d

	
#test_proteinatlas_db_dict_gen_PROTEIN()