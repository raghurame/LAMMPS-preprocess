import time as t
import decimal
import math as m
from subprocess import call

"""

REQUIRED FUNCTIONS:
~~~~~~~~~~~~~~~~~~

1. Split data file based on atom index
2. Split data file based on X, Y, Z bounds
3. Merge and dock a single chain over a substrate (input: user specifies side of substrate where chain can be docked)

"""

class lammpsdata(object):
	"""docstring for lammpsdata"""
	def __init__(self, filename):

		# from load ()
		self.filename = filename
		self.massCheck = 0
		self.atomCheck = 0
		self.bondCheck = 0
		self.angleCheck = 0
		self.dihedralCheck = 0
		self.improperCheck = 0
		self.natoms = 0
		self.nbonds = 0
		self.nangles = 0
		self.ndihedrals = 0
		self.nimpropers = 0
		self.atomTypes = 0
		self.bondTypes = 0
		self.angleTypes = 0
		self.dihedralTypes = 0
		self.improperType = 0
		self.xlo = 0
		self.xhi = 0
		self.ylo = 0
		self.yhi = 0
		self.zlo = 0
		self.zhi = 0
		self.masses = {}
		self.atom = {}
		self.bond = {}
		self.angle = {}
		self.dihedral = {}

		# from bounds ()
		self.xloBounds = 0
		self.xhiBounds = 0
		self.yloBounds = 0
		self.yhiBounds = 0
		self.zloBounds = 0
		self.zhiBounds = 0

		# from com ()
		self.xCOM = 0
		self.yCOM = 0
		self.zCOM = 0

	def load (self):

		def extract_numbers (line):
			lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
			for item in lineString.split (","):
				try:
					yield decimal.Decimal (item)
				except:
					pass

		with open (self.filename, "r") as file:

			for line in file:

				if (line.count ("atoms") > 0):
					atomInfo = list (extract_numbers (line))
					self.natoms = int (atomInfo[0])

				if (line.count ("bonds") > 0):
					bondInfo = list (extract_numbers (line))
					self.nbonds = int (bondInfo[0])

				if (line.count ("angles") > 0):
					angleInfo = list (extract_numbers (line))
					self.nangles = int (angleInfo[0])

				if (line.count ("dihedrals") > 0):
					dihedralInfo = list (extract_numbers (line))
					self.ndihedrals = int (dihedralInfo[0])

				if (line.count ("impropers") > 0):
					improperInfo = list (extract_numbers (line))
					self.nimpropers = int (improperInfo[0])

				if (line.count ("atom types") > 0):
					atomType = list (extract_numbers (line))
					self.atomTypes = int (atomType[0])

				if (line.count ("bond types") > 0):
					bondType = list (extract_numbers (line))
					self.bondTypes = int (bondType[0])

				if (line.count ("angle types") > 0):
					angleType = list (extract_numbers (line))
					self.angleTypes = int (angleType[0])

				if (line.count ("dihedral types") > 0):
					dihedralType = list (extract_numbers (line))
					self.dihedralTypes = int (dihedralType[0])

				if (line.count ("improper types") > 0):
					improperType = list (extract_numbers (line))
					self.improperTypes = int (improperType[0])

				if (line.count ("xlo xhi") > 0):
					dimX = list (extract_numbers (line))
					try:
						self.xlo = dimX[0]
						self.xhi = dimX[1]
					except:
						self.xlo = 0
						self.xhi = 0

				if (line.count ("ylo yhi") > 0):
					dimY = list (extract_numbers (line))
					try:
						self.ylo = dimY[0]
						self.yhi = dimY[1]
					except:
						self.ylo = 0
						self.yhi = 0

				if (line.count ("zlo zhi") > 0):
					dimZ = list (extract_numbers (line))
					try:
						self.zlo = dimZ[0]
						self.zhi = dimZ[1]
					except:
						self.zlo = 0
						self.zhi = 0

				if (line.count ("Masses") > 0):
					
					self.masses = {}
					self.massCheck = 1

				if (self.massCheck):
					
					try:
						massInfo = list (extract_numbers (line))
						self.masses[int (massInfo[0])] = massInfo[1]
						if (int (massInfo[0]) == int (self.atomTypes)):
							self.massCheck = 0
					except:
						pass

				if (line.count ("Atoms") > 0):
					self.atomCheck = 1
					self.atom = {}

				if (self.atomCheck == 1):

					try:
						atomLineString = list (extract_numbers (line))
						self.atom [int (atomLineString[0])] = atomLineString
						if (int (atomLineString[0]) == self.natoms):
							self.atomCheck = 0
					except:
						pass

				if (line.count ("Bonds") > 0):
					self.bondCheck = 1
					self.bond = {}

				if (self.bondCheck == 1):

					try:
						bondLineString = list (extract_numbers (line))
						self.bond [int (bondLineString[0])] = bondLineString
						if (int (bondLineString[0]) == self.nbonds):
							self.bondCheck = 0
					except:
						pass

				if (line.count ("Angles") > 0):
					self.angleCheck = 1
					self.angle = {}

				if (self.angleCheck == 1):

					try:
						angleLineString = list (extract_numbers (line))
						self.angle [int (angleLineString[0])] = angleLineString
						if (int (angleLineString[0]) == self.nangles):
							self.angleCheck = 0
					except:
						pass

				if (line.count ("Dihedrals") > 0):
					self.dihedralCheck = 1
					self.dihedral = {}

				if (self.dihedralCheck == 1):

					try:
						dihedralLineString = list (extract_numbers (line))
						self.dihedral [int (dihedralLineString[0])] = dihedralLineString

						if (int (dihedralLineString[0]) == self.ndihedrals):
							self.dihedralCheck = 0
					except:
						pass


					# break

				# if (line.count ("Masses") > 0):
				# 	improperInfo = list (extract_numbers (line))
				# 	self.nimpropers = int (improperInfo[0])

				# while (line.count ("Atoms") > 0):
				# 	pass
				# dataList = list (extract_numbers (line))
				# # print (len (dataList))
				# # print (dataList[0]+3)
				# t.sleep(1)
		
	def bounds (self):

		for atomID, info in self.atom.items ():

			if (int (info[3]) < self.xloBounds):
				self.xloBounds = info[3]
			if (int (info[3]) > self.xhiBounds):
				self.xhiBounds = info[3]

			if (int (info[4]) < self.yloBounds):
				self.yloBounds = info[4]
			if (int (info[4]) > self.yhiBounds):
				self.yhiBounds = info[4]

			if (int (info[5]) < self.zloBounds):
				self.zloBounds = info[5]
			if (int (info[5]) > self.zhiBounds):
				self.zhiBounds = info[5]

	def com (self):

		for atomID, info in self.atom.items ():

			self.xCOM += info[3]
			self.yCOM += info[4]
			self.zCOM += info[5]

		self.xCOM /= self.natoms
		self.yCOM /= self.natoms
		self.zCOM /= self.natoms

	def recenter(self, dx = 0, dy = 0, dz = 0):
		
		for atomID, info in self.atom.items ():

			info[3] += dx
			info[4] += dy
			info[5] += dz

	def split (self, )

class betaipp(object):
	"""docstring for betaipp"""
	def __init__(self, filename):
		self.filename = filename
		self.removeID = []
		self.cleaned = {}

	def clean (self, deleteCoords):

		betaipp_data = lammpsdata (self.filename)
		betaipp_data.load ()

		deleteX = []
		deleteY = []
		newID = 1

		with open (deleteCoords, "r") as file:

			for line in file:

				if (line.count ("Info) x: ") > 0):
					
					X = (line.replace ("Info) x: ", "").replace ("\n", "").strip ())
					deleteX.append ( float (X))

				if (line.count ("Info) y: ") > 0):
					
					Y = (line.replace ("Info) y: ", "").replace ("\n", "").strip ())
					deleteY.append (float (Y))

		for atomID, info in betaipp_data.atom.items ():

			for i in range (0, len (deleteX)):

				distance = m.sqrt (m.pow ((float (info[3]) - deleteX[i]), 2) + m.pow ((float (info[4]) - deleteY[i]), 2))

				if (distance < 1):
					
					self.removeID.append (atomID)

		for atomID, info in betaipp_data.atom.items ():

			if atomID not in self.removeID:

				self.cleaned [newID] = info
				newID += 1

		with open ("cleanedBeta.data", "w") as file:
			
			for atomID, data in self.cleaned.items ():

				file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (atomID, data[1], data[2], round (data[3], 4), round (data[4], 4), round (data[5], 4)))

	def recenter (self, outputFile = "recenteredBeta.data"):

		beta = lammpsdata (self.filename)
		beta.load ()
		beta.com ()
		beta.recenter (dx = -beta.xCOM, dy = -beta.yCOM, dz = -beta.zCOM)
		beta.com ()
		beta.bounds ()

		# printing the recentered data into new file
		with open (outputFile, "w") as file:

			file.write ("Created by preprocess module in python\n\n")
			file.write ("\t{}\tatoms\n".format (beta.natoms))
			file.write ("\t{}\tbonds\n".format (beta.nbonds))
			file.write ("\t{}\tangles\n".format (beta.nangles))
			file.write ("\t{}\tdihedrals\n".format (beta.ndihedrals))
			file.write ("\t{}\timpropers\n".format (beta.nimpropers))

			file.write ("\n")

			file.write ("\t{}\tatom types\n".format (beta.atomTypes))
			file.write ("\t{}\tbond types\n".format (beta.bondTypes))
			file.write ("\t{}\tangle types\n".format (beta.angleTypes))
			file.write ("\t{}\tdihedral types\n".format (beta.dihedralTypes))
			file.write ("\t{}\timproper types\n".format (beta.improperType))

			file.write ("\n")

			file.write ("\t{}\t{}\txlo\txhi\n".format (round (beta.xloBounds, 4), round (beta.xhiBounds, 4)))
			file.write ("\t{}\t{}\tylo\tyhi\n".format (round (beta.yloBounds, 4), round (beta.yhiBounds, 4)))
			file.write ("\t{}\t{}\tzlo\tzhi\n".format (round (beta.zloBounds, 4), round (beta.zhiBounds, 4)))

			file.write ("\n")

			file.write ("Masses\n\n")

			for x in range (1, beta.atomTypes+1):

				file.write ("\t{}\t{}\n".format (x, beta.masses[x]))

			file.write ("\nAtoms\n\n")

			for atomID, data in beta.atom.items ():

				file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (data[0], data[1], data[2], round (data[3], 4), round (data[4], 4), round (data[5], 4)))

		call (["perl", "lammps2pdb.pl", outputFile.replace (".data", "")])

