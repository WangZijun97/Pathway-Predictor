import requests
import urllib.parse
from functools import reduce
from settings import THRESHOLD
from datetime import datetime

class TissueType:
	def __init__(self, name, colcode):
		self.name = name
		self.col = colcode
		
tissuetypes = [
	"adipose tissue",
	"adrenal gland",
	"amygdala",
	"appendix",
	"basal ganglia",
	"bone marrow",
	"breast",
	"cerebellum",
	"cerebral cortex",
	#cervix, uterine removed due to non-functional url link
	#"cervix, uterine",
	"colon",
	"corpus callosum",
	"ductus deferens",
	"duodenum",
	"endometrium 1",
	"epididymis",
	"esophagus",
	"esophagus",
	"fallopian tube",
	"gallbladder",
	"heart muscle",
	"hippocampal formation",
	"hypothalamus",
	"kidney",
	"liver",
	"lung",
	"lymph node",
	"midbrain",
	"olfactory region",
	"ovary",
	"pancreas",
	"parathyroid gland",
	"placenta",
	"pons and medulla",
	"prostate",
	"rectum",
	"retina",
	"salivary gland",
	"seminal vesicle",
	"skeletal muscle",
	"skin 1",
	"small intestine",
	"smooth muscle",
	"spinal cord",
	"spleen",
	"stomach 1",
	"testis",
	"thalamus",
	"thymus",
	"thyroid gland",
	"tongue",
	"tonsil",
	"urinary bladder",
	"vagina",
	"B-cells",
	"dendritic cells",
	"granulocytes",
	"monocytes",
	"NK-cells",
	"T-cells",
	"total PBMC"
	]
	
#convert raw str type tissue name to param for proteinatlas api url	
def encodetissue(tissue):
	return "t_RNA_" + tissue.replace(" ", "_")
	
def decodecol(colheader):
	return colheader[13:-5]
	
tissuetypesasparam = reduce(lambda x, y: x + "," + y,
	map(lambda x: encodetissue(x), tissuetypes))

#testing validity of tissuetypes	
def tissuetypestest(list, encode = True):
	for i in list:
		if encode:
			code = encodetissue(i)
		else:
			code = i
		#url = "https://www.proteinatlas.org/api/search_download.php?search=ENSG00000178796&format=json&columns={0}&compress=no".format(code)
		
		params = {"search": "ENSG00000178796", "format": "json", "columns": code, "compress": "no"}
		urlparams = urllib.parse.urlencode(params, doseq=True)
		url = "https://www.proteinatlas.org/api/search_download.php?" + urlparams
		
		data = requests.get(url).json()
		print("{0} > {2} | {1}".format(i, data, url))

#given a protein code, return set of tissue types in raw str format		
def gettissues(_logger, protein, mem, threshold = THRESHOLD):
	
	#check memoi table
	if protein in mem:
		return set(mem[protein])
		
	else:
		print(_logger.log("[{0}] Getting Tissue Types of {1} with NX >= {2}".format(datetime.now(), protein, threshold)))
		#convert ensp to ensg code
		#code = ensmblptog(protein)
		code = protein
		tissueset = set()
		
		#request proteinatlas api for data
		params = {"search": code, "format": "json", "columns": tissuetypesasparam, "compress": "no"}
		urlparams = urllib.parse.urlencode(params, doseq=True)
		url = "https://www.proteinatlas.org/api/search_download.php?" + urlparams
		data = requests.get(url).json()
		if len(data) < 1:
			print(_logger.log("{0} | {1}".format(code, data)))
			mem[protein] = []
			return tissueset
		else:
			data = data[0]
		
		#testing
		"""
		for i in data:
			print(" > {0}: {1}".format(decodecol(i), data[i]))
			"""
		
		#filter all NX < 1.0
		for i in data:
			try:
				if float(data[i]) >= threshold:
					tissueset.add(decodecol(i))
			except:
				continue
		
		mem[protein] = list(tissueset)
		return tissueset

#checks if code entered is Ensembl Gene ID
#used in ensmblptog()	
def isgene(code):
	if len(code) < 4:
		return False
	if code[0:4] == "ENSG":
		return True
	return False

#converts ENSP to ENSG. Currently unused	
def ensmblptog(_logger, proteincode):
	print(_logger.log("[{0}] Converting ENSP to ENSG for {1}".format(datetime.now(), proteincode)))
	
	def urlhelper(code):
		return "http://rest.ensembl.org/lookup/id/{0}.json".format(code)
		
	data = requests.get(urlhelper(proteincode)).json()
	#print(proteincode, data)
	temp = data["Parent"]

	while not isgene(temp):
		data = requests.get(urlhelper(temp)).json()
		#print(temp, data)
		temp = data["Parent"]
		
	print(" > " + temp)
	return temp
		
#tissuetypestest(tissuetypes, encode = True) #passed
#print(ensmblptog("ENSP00000373702")) #passed
#print(ensmblptog("ENSP00000216274")) #passed
#print(ensmblptog("ENSP00000279441")) #passed
#print(tissuetypesasparam) #looks right
#tissuetypestest([tissuetypesasparam], encode = False) #passed
#print(gettissues("ENSP00000216274", {})) #passed
		
		