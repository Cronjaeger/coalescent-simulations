#### NO LONGER IN USE
#### ALL CODE IN ilbCoal.py

class PPP(object):
	'A Possion Point Process, for driving Xi-coalescents'
	
	def __init__(self,intensityMeasure,n,T):
		self.measure = intensityMeasure
		self.events = {}
		self.resample()
		
	#def resample(self):
		#while t <= T:
			#do stuff
	
	def __getitem__(self):
		#do stuff

	#def TMRCA(self):
		#if self.hasCoalesced:
			#return max(self.events.keys())
		#else return nan
