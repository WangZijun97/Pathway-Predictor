import datetime

class Logger:
	def __init__(self, time):
		self.logls = []
		self.name = str(time).replace(":", ".")
		
	def log(self, str, toprint = False, timestamp = False):
		if timestamp:
			str = "[{0}] ".format(datetime.datetime.now()) + str
		self.logls.append(str)
		if toprint:
			print(str)
		return str
		
	def printout(self):
		return "\n".join(self.logls)