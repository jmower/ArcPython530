import numpy as np
class Vertex:
	def __init__(self, easting, elevation, northing):
		self.easting = easting
		self.elevation = elevation
		self.northing = northing
		self.reserve = 1
		self.accum = 0

# v = Vertex(580000,100,4720000)
# s = 'E ' + repr(v.easting) + ' N ' + repr(v.northing ) + ' Elev ' + repr(v.elevation) + ' Reserve ' + repr(v.reserve) + ' Accum ' + repr(v.accum)
# print s

