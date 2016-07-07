import Vertex
import re # regular expressions
import numpy as np
class ArbGrid:
	# constants
	NURBSDEGREE = 3
	P = NURBSDEGREE
	ORDER = P + 1
	BASISVALCUTOFF = 0.0000001
	EVALSAMPLINGRESDEFAULT = 100 # 25 meter sampling res in u and v for evaluating points to surface
	tessVertices = [] # list of evaluated vertices arranged by columns. Each column is a list of Vertex objects. Columns increase in easting by
						# the given evaluated sampling distance. Within a column, samples increase with northing by the evaluated sampling distance
	tessTriangles = [] # list of evaluated triangles
	tessTriNormals = []

	def __init__(self, agFilename):
		self.agFile = open(agFilename, 'r')
		self.readHeader()
		self.readData()
		self.buildBSpline()

	def readHeader(self):
		# ArbGrid reading constants
		MINEASTING_START = 0
		MINNORTHING_START = 9
		MAXEASTING_START = 18
		MAXNORTHING_START = 27
		ROWSAMPDIST_START = 36
		COLSAMPDIST_START = 45
		NUMCOLS_START = 54
		NUMROWS_START = 59
		DATA_START = 64
		ELEV_LENGTH = 6
		# grab the entire ArbGrid
		self.header = self.agFile.readline()
		# get substrings for header information and for the remainder of the elevation data
		sMinEasting = self.header[MINEASTING_START:MINNORTHING_START]
		sMinNorthing = self.header[MINNORTHING_START:MAXEASTING_START]
		sMaxEasting = self.header[MAXEASTING_START:MAXNORTHING_START]
		sMaxNorthing = self.header[MAXNORTHING_START:ROWSAMPDIST_START]	
		sRowSampDist = self.header[ROWSAMPDIST_START:COLSAMPDIST_START]
		sColSampDist = self.header[COLSAMPDIST_START:NUMCOLS_START]
		sNumCols = self.header[NUMCOLS_START:NUMROWS_START]
		sNumRows = self.header[NUMROWS_START:DATA_START]
		self.sElevations = self.header[DATA_START:]
		self.minEasting = float(sMinEasting)
		self.minNorthing = float(sMinNorthing)
		self.maxEasting = float(sMaxEasting)
		self.maxNorthing = float(sMaxNorthing)
		self.rowRes = float(sRowSampDist)
		self.colRes = float(sColSampDist)
		self.numRows = int(sNumRows)
		self.numCols = int(sNumCols)

		# List comprehension for extracting floating point elevations from each of the string versions (self.sElevations) of the elevation data elements
		self.elevations = [float(self.sElevations[i * ELEV_LENGTH : i * ELEV_LENGTH + ELEV_LENGTH]) for i in range(0, self.numRows * self.numCols)]
		# print "self.elevations:", self.elevations
	def readData(self):
		# List comprehension. For each easting, elevation, northing tuple, write the tuple to a list as a Vertex object
		# backslash is inewrap character
		# self.vertices = [Vertex.Vertex(self.minEasting + easting * self.colRes, self.elevations[northing * self.numCols + easting], \
		# 	self.minNorthing + northing * self.rowRes) for easting in range(0, self.numCols) for northing in range(0, self.numRows)]
		self.vertices = [Vertex.Vertex(self.minEasting + easting * self.colRes, self.elevations[easting * self.numRows + northing], \
			self.minNorthing + northing * self.rowRes) for easting in range(0, self.numCols) for northing in range(0, self.numRows)]
		
		# print "readData, values in self.vertices:" # prints all northings in an easting, followed by next easting but elevations are transposed by row/col.
		# for vert in self.vertices:
			# print (vert.easting, vert.elevation, vert.northing)

	def buildBSpline(self):
		# most of the code for building the b-spline comes from PenAndInkNURBS
		# print self.NURBSDEGREE
		if (self.numCols < self.NURBSDEGREE or self.numRows < self.NURBSDEGREE):
			return False
		self.numberOfArbGridUKnots =self.numCols + self.P + 1
		self.numberOfArbGridVKnots = self.numRows + self.P + 1
		self.uKnots = []
		self.vKnots = []
		self.buildUniformArbGridKnotArrayForArbGrid()

	def buildUniformArbGridKnotArrayForArbGrid(self):
		# find the knot sequence in the U direction
		iM = self.numberOfArbGridUKnots;
		lowRange = self.ORDER
		midRange = iM - self.ORDER
		hiRange = iM
		for k in range(0, lowRange):
			self.uKnots.append(0.0)
		for k in range(lowRange, midRange):
			self.uKnots.append(float(k - (self.ORDER - 2)))
		for k in range(midRange, hiRange):
			self.uKnots.append(float(iM - (self.ORDER + 1)))
		# print "buildUniformArbGridKnotArrayForArbGrid: uKnots", self.uKnots
		# now find the knot sequence in the V direction
		iM = self.numberOfArbGridVKnots;
		midRange = iM - self.ORDER
		hiRange = iM
		for k in range(0, lowRange):
			self.vKnots.append(0.0)
		for k in range(lowRange, midRange):
			self.vKnots.append(float(k - (self.ORDER - 2)))
		for k in range(midRange, hiRange):
			self.vKnots.append(float(iM - (self.ORDER + 1)))
		# print "buildUniformArbGridKnotArrayForArbGrid: vKnots", self.vKnots

	def evalVertexOnNurbsSurface(self, easting, northing):
		# print "evalVertexOnNurbsSurface: evaluating" + repr(easting) + ", " + repr(northing)
		minU = self.uKnots[self.NURBSDEGREE - 1]
		minV = self.vKnots[self.NURBSDEGREE - 1]
		#print minU, minV
		maxU = self.uKnots[len(self.uKnots) - (self.NURBSDEGREE - 1)]
		maxV = self.vKnots[len(self.vKnots) - (self.NURBSDEGREE - 1)]
		#print maxU, maxV
		u = (maxU - minU) * (easting - self.minEasting) / (self.maxEasting - self.minEasting);
		v = (maxV - minV) * (northing - self.minNorthing) / (self.maxNorthing - self.minNorthing);
		# print "evalVertexOnNurbsSurface:", u, v, easting, northing, self.minEasting, self.minNorthing, self.maxEasting, self.maxNorthing
		#print u, v
		vertex = self.evalVertex(u, v)
		# print "evalVertexOnNurbsSurface: evaluated vertex: ", vertex.easting, vertex.northing
		return vertex

	def evalVertex(self, u, v):
		colVect = []
		for row in range(0, self.numRows):
			colVect.append(self.evalVertexU(u, row)) # colVect values look the same as for PenAndInkNURBS
		# print "evalVertex: colVect:"
		# for aVert in colVect:
		# 	print "evalVertex: aVert:", aVert.easting, aVert.elevation, aVert.northing # ok
		vertex = self.evalVertexV(v, colVect)
		# print "evalVertex: evaluated vertex: " + repr(vertex) # ok
		return vertex

	def evalVertexU(self, u, row):
		
		vertex = Vertex.Vertex(0, 0, 0)
		for col in range(0, self.numCols):
			basisVal = self.CoxDeBoor(u, col, self.NURBSDEGREE + 1, self.uKnots)
			# print "evalVertexU, basisVal:", basisVal
			cpIndex = col * self.numRows + row
			if basisVal > self.BASISVALCUTOFF:
				# np.array([1, 2, 3])
				# print "unevaluated vertex in evalVertexU:", self.vertices[cpIndex].easting, self.vertices[cpIndex].elevation, self.vertices[cpIndex].northing # correct values printing
				easting = basisVal * self.vertices[cpIndex].easting
				elevation = basisVal * self.vertices[cpIndex].elevation
				northing = basisVal * self.vertices[cpIndex].northing
				#vertex = np.array([easting, elevation, northing])
				vertex.easting += easting
				vertex.elevation += elevation
				vertex.northing += northing
				#print "evaluated vertex in evalVertexU:", vertex.easting, vertex.elevation, vertex.northing, basisVal # reasonable output
		# print "evalVertexU: returning vertex: ", vertex.easting, vertex.elevation, vertex.northing # ok
		return vertex

	def evalVertexV(self, v, colVect):
		vertex = Vertex.Vertex(0, 0, 0)
		# for aVert in colVect:
		# 	print "evalVertexV: aVert:", aVert.easting, aVert.elevation, aVert.northing # ok
		for row in range(0, self.numRows):
			basisVal = self.CoxDeBoor(v, row, self.NURBSDEGREE + 1, self.vKnots)
			if basisVal > self.BASISVALCUTOFF:
				easting = basisVal * colVect[row].easting
				elevation = basisVal * colVect[row].elevation
				northing = basisVal * colVect[row].northing
				# print "evalVertexV:", easting, elevation, northing
				#vertex = Vertex.Vertex(easting, elevation, northing)
				vertex.easting += easting;
				vertex.elevation += elevation
				vertex.northing += northing
		# print "velVertexV: returning ", vertex.easting, vertex.elevation, vertex.northing
		return vertex

	def CoxDeBoor(self, u, i, k, knots):
		if (k == 1):
			if knots[i] <= u and u <= knots[i +  1]:
				return 1.0
			return 0.0
		den1 = knots[i + k - 1] - knots[i]
		den2 = knots[i + k] - knots[i + 1]
		eq1 = 0.0
		eq2 = 0.0
		if den1 > 0.0:
			eq1 = ((u - knots[i]) / den1) * self.CoxDeBoor(u, i, k-1, knots)
	 	if den2 > 0.0:
	 		eq2 = (knots[i + k] - u) / den2 * self.CoxDeBoor(u, i+1, k-1, knots);
	 	return eq1 + eq2

	def buildTessTriangles(self):
		# print "in buildTessTriangles"
		# Use the vertices in self.tessVertices to build triangles. Each triangle will be composed as a list of 3 indexes into self.tessVertices. Each
		# list of vertices will be held in self.tessTriangles.
		# Each triangle will also have a computed surface normal to be stored in self.tessTriNormals by index into self.tessTriangles.
		# Vertex indices are stored in clockwise order for each triple.
		############# U = upper tri
		#         # # L = lower tri 
		#   U   #   # Vertex 0 is at southwest corner. Vertex 1 is one sample north of 0. Vertex numRows is one vertex east of Vertex 0
		#     #     #
		#   #   L   #
		# #         #
		#############
		for col in range(0, self.numCols - 1):
			for row in range(0, self.numRows - 1):
				# compute a linear index
				index = col * self.numRows + row
				# lower triangle
				i0 = index
				i1 = index + self.numRows + 1
				i2 = index + self.numRows
				self.tessTriangles.append([i0, i1, i2])
				# upper triangle
				i0 = index
				i1 = index + 1
				i2 = index + 1 + self.numRows
				self.tessTriangles.append([i0, i1, i2])

		# now create normals
		for tri in self.tessTriangles:
			# each tri is a list of indices into self.tess.vertices
			linearIndex0 = tri[0]
			linearIndex1 = tri[1]
			linearIndex2 = tri[2]
			linearIndices = [linearIndex0, linearIndex1, linearIndex2]
			# convert raw index into one that finds the correct colVect in self.tessVerts and then the correct row in the colVect
			triCoords = []
			for li in linearIndices:				
				col = li / self.numRows
				row = li % self.numRows
				vert = self.tessVertices[col][row]
				triCoords.append(np.array([vert.easting, vert.elevation, vert.northing]))
			#print "buildTessTriangles: triCoords:", triCoords
			v0 = triCoords[0] - triCoords[1]
			v1 = triCoords[1] - triCoords[2]
			normal = np.cross(v0,v1)
			magnitude = np.linalg.norm(normal)
			normal /= magnitude # normalize vector components
			self.tessTriNormals.append(normal)
		# for norm in self.tessTriNormals:
		# 	print "builTessTriangles, norm:", norm
		#print "buildTessTriangles: normal: ", normal
	def evaluateSurface(self):
		for easting in range (int(self.minEasting), int(self.maxEasting), self.EVALSAMPLINGRESDEFAULT):
			colVect = []
			for northing in range(int(self.minNorthing), int(self.maxNorthing), self.EVALSAMPLINGRESDEFAULT):
				colVect.append(self.evalVertexOnNurbsSurface(easting, northing))
			self.tessVertices.append(colVect)
		
		self.buildTessTriangles()

		# for tri in self.tessTriangles: # just for debugging
		# 	print tri 
		# for col in self.tessVertices: # just for debugging
		# 	for vert in col:
		# 		print vert.easting, vert.elevation, vert.northing