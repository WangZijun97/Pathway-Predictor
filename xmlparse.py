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
	