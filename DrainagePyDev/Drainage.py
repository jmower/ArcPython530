import ArbGrid as ag
import numpy as np
from Tkinter import *
from tkFileDialog import askopenfilename

# this code was modified from http://effbot.org/tkinterbook/tkinter-hello-again.htm
class App:
	sunVector = np.array([[-.707106], [.707106], [0]])
	farZ = 10000.0
	nearZ = 100.0
	theArbrid = None
	drawnObjects = []
	# v = None #Experiment--global variables for writable label (through v as a StringVar())
	# textBox = None
	def __init__(self, master):
		#print "in init"
		frame = Frame(master)
		frame.pack()

		self.quitButton = Button(
		    frame, text="QUIT", fg="red", command=frame.quit
		    )
		self.quitButton.pack(side=LEFT)

		self.LoadAgButton = Button(frame, text="Build ArbGrid", command=lambda : self.buildArbGrid()) # if I don't make this a lambda function
																									# it buildArbGrid() executes without button push
		self.LoadAgButton.pack(side=LEFT)

		entryOptions = {} # create a dictionary and add options to it below. '**' notation below allows 'kwargs' (keyword args)
		entryOptions['width'] = '10'

		label = Label(frame, text="VP E")
		label.pack(side=LEFT)
		self.vpEastingEntry = Entry(frame, entryOptions)
		self.vpEastingEntry.insert(0, "580400")
		self.vpEastingEntry.pack(side=LEFT)

		label = Label(frame, text="VP N")
		label.pack(side=LEFT)
		self.vpNorthingEntry = Entry(frame, entryOptions)
		self.vpNorthingEntry.insert(0, "4717000")
		self.vpNorthingEntry.pack(side=LEFT)

		label = Label(frame, text="VP El")
		label.pack(side=LEFT)
		entryOptions['width'] = '5'
		self.vpElevationEntry = Entry(frame, entryOptions)
		self.vpElevationEntry.insert(0, "1000")
		self.vpElevationEntry.pack(side=LEFT)

		label = Label(frame, text="VP Az")
		label.pack(side=LEFT)
		self.vpAzimuthEntry = Entry(frame, entryOptions)
		self.vpAzimuthEntry.insert(0, "0")
		self.vpAzimuthEntry.pack(side=LEFT)

		label = Label(frame, text="VP Alt")
		label.pack(side=LEFT)
		self.vpAltitudeEntry = Entry(frame, entryOptions)
		self.vpAltitudeEntry.insert(0, "10")
		self.vpAltitudeEntry.pack(side=LEFT)

		label = Label(frame, text="VP fovy")
		label.pack(side=LEFT)
		self.fovyEntry = Entry(frame, entryOptions)
		self.fovyEntry.insert(0, "22.0")
		self.fovyEntry.pack(side=LEFT)

		# Experiment--make a writable label--works
		# self.v = StringVar()
		# self.textBox = Label(frame, textvariable=self.v)
		# self.v.set("more junk")
		# self.textBox.pack(side=LEFT)

		# Create a canvas widget
		self.canvas = Canvas(master, width=512, height=256)
		self.canvas.pack()
	def loadOptionsFromEntries(self):
		# parse text entry values. Need to wrap in try...except blocks
		vpEastingS = self.vpEastingEntry.get()
		vpEasting = float(vpEastingS)
		vpNorthingS = self.vpNorthingEntry.get()
		vpNorthing = float(vpNorthingS)
		vpElevationS = self.vpElevationEntry.get()
		vpElevation = float(vpElevationS)
		vpAzimuthS = self.vpAzimuthEntry.get()
		vpAzimuthDD = float(vpAzimuthS)
		vpAzimuth = np.radians(vpAzimuthDD)
		vpAltitudeS = self.vpAltitudeEntry.get()
		vpAltitudeDD = float(vpAltitudeS)
		vpAltitude = np.radians(vpAltitudeDD)
		fovyS = self.fovyEntry.get()
		fovyDD = float(fovyS)
		fovy = np.radians(fovyDD)
		self.viewpoint = np.array([[vpEasting], [vpElevation], [vpNorthing], [1.0]])
		self.azimuth = vpAzimuth
		self.altitude = vpAltitude
		self.fovy = fovy
		print vpEasting
	def buildTransforms(self):
		# parse text entry values
		self.loadOptionsFromEntries()

		windowWidth = self.canvas.winfo_width()
		windowHeight = self.canvas.winfo_height()
		tX = -self.viewpoint[0][0]
		tY = -self.viewpoint[1][0]
		tZ = -self.viewpoint[2][0]
		sX = 1
		sY = -1 # flip y axis
		self.model = np.array([[1,0,0,tX],[0,1,0,tY],[0,0,1,tZ],[0,0,0,1]])

		# adapt code from C:\Users\jmower\Documents\Urhere\src\ParallelDB\PenAndInkNURBS\OpenGLControl.cpp::buildTransforms()
		# will need to find equivalent of glm::perspective()
		rotAzimuthDD = np.degrees(self.azimuth) + 180.0
		if (rotAzimuthDD > 360.0):
			rotAzimuthDD = rotAzimuthDD % 360.0
		rotAzimuth = rotAzimuthDD / 180.0 * np.pi

		# follows ParallelDrainage.COpenGLControl::buildTransforms() in 'test' section. Using column major format here (same as used by glm and OpenGL)
		identity = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
		self.view = identity

		scaleM = np.array([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
		self.view = np.dot(self.view, scaleM)

		rotAltitudeM = np.array([[1,0,0,0],[0,np.cos(self.altitude),-np.sin(self.altitude),0],[0,np.sin(self.altitude),np.cos(self.altitude),0],[0,0,0,1]]) # rotate about X axis
		self.view = np.dot(self.view, rotAltitudeM)

		rotAzimuthM = np.array([[np.cos(-rotAzimuth),0,np.sin(-rotAzimuth),0],[0,1,0,0],[-np.sin(-rotAzimuth),0,np.cos(-rotAzimuth),0],[0,0,0,1]]) # rotate about Y axis
		self.view = np.dot(self.view, rotAzimuthM)

		aspect = float(self.canvas.winfo_width()) / float(self.canvas.winfo_height())

		modelView = np.dot(self.view, self.model)

		f = 1.0/np.tan(self.fovy/2)

		self.perspective = np.array([[f/aspect,0,0,0],[0,f,0,0],[0,0,(self.farZ+self.nearZ)/(self.nearZ-self.farZ),(2*self.farZ*self.nearZ)/(self.nearZ-self.farZ)],[0,0,-1,0]])

		self.mvp = np.dot(self.perspective, modelView)

	def drawSurface(self):

		#implement Z buffering??

		# get local variable names for lists in theArbGrid
		tessT = self.theArbGrid.tessTriangles # list of triangles with 3 vertex indices (into tessVertices) per triangle record
		tessV = self.theArbGrid.tessVertices # list of all tessellated vertices, stored as on Vertex object per record
		tessN = self.theArbGrid.tessTriNormals # surface normals for each triangle in tessTriangles

		t = 0
		for tri in tessT:
			# get linear indices into tessV for one triangle
			linearIndex0 = tri[0]
			linearIndex1 = tri[1]
			linearIndex2 = tri[2]
			linearIndices = [linearIndex0, linearIndex1, linearIndex2]
			# convert raw index into one that finds the correct colVect in self.tessVerts and then the correct row in the colVect
			triCoords = []
			for li in linearIndices:	# fill triCoords with the 3 Vertex objects pointed to by linearIndices
				col = li / self.theArbGrid.numRows
				row = li % self.theArbGrid.numRows
				vert = self.theArbGrid.tessVertices[col][row]
				# print vert.easting, vert.elevation, vert.northing
				triCoords.append(vert)
			triNorm = tessN[t]

			t += 1

			# find the angle between this triangle's normal and the incoming sun vector.
			sunVTranspose = np.transpose(self.sunVector)
			cosNormSun = np.dot(sunVTranspose, triNorm)

			# Attenuate surface color by cosNormSun
			colorRGB = np.array([0, 255, 0])
			attenuatedColorRGB = colorRGB * np.abs(cosNormSun) # don't want negative values. Probably want to have minimum attenuated color shading (> 0)
			attenuatedColorHex = self.makeHexColorString(attenuatedColorRGB)

			# project vertices by transformation matrix
			tct = []
			windowWidth = self.canvas.winfo_width()
			windowHeight = self.canvas.winfo_height()
			for vert in triCoords:
				# project each vertex in triCoords my mvp
				v = np.array([[vert.easting], [vert.elevation], [vert.northing], [1]])
				vertTrans = np.dot(self.mvp, v) # transform vertTrans to clip coords

				vertTrans = vertTrans / vertTrans[3] # transform vertTrans to normalized device coords

				x = y = 0.0 # coordinates of viewport
				# transform vertTrans x and y to window coords: http://home.deec.uc.pt/~peixoto/eda/opengl/glViewport.html. Z coord holds depth (vertTrans[2][0]) and should not be scaled.
				vertTrans[0][0] = (vertTrans[0][0] + 1) * (windowWidth / 2.0) + x
				vertTrans[1][0] = windowHeight - ((vertTrans[1][0] + 1) * (windowHeight / 2.0) + y)
				tct.append(vertTrans)

			# draw the 3 projected vertices for this triangle as a triangular polygon
			self.drawnObjects.append(self.canvas.create_polygon(int(tct[0][0]), int(tct[0][1]), int(tct[1][0]), int(tct[1][1]), int(tct[2][0]), int(tct[2][1]), fill=attenuatedColorHex, outline='yellow'))

		# experiment--write to text label
		#self.v = StringVar()
		# self.v.set("finished processing")
	def makeHexColorString(self, attenuatedColorRGB):
		r, g, b = attenuatedColorRGB
		hexR = hex(int(r))[2:] # cut off leading 0x
		hexG = hex(int(g))[2:]
		hexB = hex(int(b))[2:]
		# print hexR, hexG, hexB
		if len(hexR) < 2:
			hexR = "0" + hexR
		if len(hexG) < 2:
			hexG = "0" + hexG
		if len(hexB) < 2:
			hexB = "0" + hexB
		colorString = "#" + hexR + hexG + hexB

		return colorString

	def buildArbGrid(self):
		# search for file
		# The following options dictionary code was modified from
		#   https://tkinter.unpythonic.net/wiki/tkFileDialog#CA-87fa57b85635fb29da2b08946d341845aecff099_19
		self.file_opt = options = {} # create a dictionary and add options to it below. '**' notation below allows 'kwargs' (keyword args)
										# of variable length (like self.file_opt) to be passed to a function (like askopenfilename() )
		options['defaultextension'] = '.ag'
		options['filetypes'] = [('all files', '.*'), ('ArbGrid files', '.ag')]
		#options['initialdir'] = '~jimmower/Documents/dev/Python/Drainage' # mac starting folder
		options['initialdir'] = 'R:\Python\Drainage' # work desktop starting folder
		options['initialfile'] = 'Angle9x9.ag'
		options['parent'] = root
		options['title'] = 'ArbGrid Files'

		filename = askopenfilename(**self.file_opt)
		if not filename:
			return

		[self.canvas.delete(id) for id in self.drawnObjects] # delete all the polygons from the previous rendering
		self.theArbGrid = ag.ArbGrid(filename)

		self.theArbGrid.evaluateSurface() # create a tessellated surface by sampling the B-spline created in theArbGrid
		self.buildTransforms() # do all the model, view, perspective, and viewport transformations
		self.drawSurface() # draw the surface

root = Tk()

app = App(root)

root.mainloop()
root.destroy() # optional; see description below
