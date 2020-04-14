import datetime
THRESHOLD = 14.8

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
	
def encodecoltissue(tissue):
	return "Tissue RNA - " + tissue + " [NX]"
	
def decodecol(colheader):
	return colheader[13:-5]
	
tissuetypesheaders = list(map(lambda x: encodecoltissue(x), tissuetypes))


#using proteinatlas tissue RNA implementation	
def gettissues(_logger, protein, db, threshold = THRESHOLD): #uniprot!!!
	#_logger.log("[{0}] Getting Tissue Types of {1} with NX >= {2}".format(datetime.datetime.now(), protein, threshold))
	
	if protein not in db:
		_logger.log(" > No Entry Found for {0}".format(protein), toprint=False)
		return set()
		
	data = db[protein]
	tissueset = set()
	
	for col in tissuetypesheaders:
		if not data[col]:
			continue
		if float(data[col]) >= threshold:
			tissueset.add(decodecol(col))
			
	return tissueset
	
	
	
#using proteinatlas protein expression implementation	
def gettissues_protein(_logger, protein, db): #ENSMBL_ID fml!!!
	#_logger.log("[{0}] Getting Tissue Types of {1} with NX >= {2}".format(datetime.datetime.now(), protein, threshold))
	
	if protein not in db:
		_logger.log(" > No Entry Found for {0}".format(protein), toprint=False)
		return set()
		
	return set(db[protein])
	
	
	
	