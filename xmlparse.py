import xml.etree.ElementTree as ElementTree
import requests

def getName(acc, mem):
	if acc in mem:
		return mem[acc]
	else:
		xml = ElementTree.fromstring(requests.get("https://www.uniprot.org/uniprot/{0}.xml".format(acc)).text)
		names = xml.iter("{http://uniprot.org/uniprot}fullName")
		name = next(names).text
		mem[acc] = name
	return name
	
def getNameFromSTRING(id, mem, mem2):
	#goes via ACC, not very reliable due to lack of perfect STRING -> UNIPROT mapping
	id = "9606"+id
	if id in mem:
		return mem[id]
	else:
		acc = requests.get("https://uniprot.org/uploadlists/?from=STRING_ID&to=ACC&format=tab&query={0}".format(id)).text.split('\n')[1].split('\t')[1]
		mem[id] = acc
		name = getName(acc, mem2)
		return name
		
def getNameFromSTRING(id, mem):
	if id in mem:
		return mem[id]
	else:
		name = requests.get("https://string-db.org/api/json/get_string_ids?identifiers={0}".format(id)).json()[0]["preferredName"]
		mem[id] = name
		return name