import numpy as np
from scipy.spatial.distance import euclidean
from math import cos, sin, pi


class Molecule:
	'''
	Class representation of a molecule
	'''

	def __init__(self, position=(0.,0.,0.), rotation=[0.,0.,0.], file=None):
		#check if file ends with xyz and try to load it
		if file is not None and file.endswith('.xyz'):
			self._load_xyz(file)

		self._bond_len_thresh = { #bond length thresholds between different atoms
			'CC': 1.60,	#atoms must be sorted alphabetically
			'CH': 1.2,
			'CO': 1.50,
			'HO': 1.00,
			'HH': 0.74,
			'OO': 1.49,
			'NO': 1.50,
			'CN': 1.49,
			'NN': 1.47,
			'HN': 1.03,
			'CFe': 0.92,
			'FeH': 0.70,
			'FeFe': 2.0,
			'MgN': 2.15,
			'CMg': 1.87,
			'HMg': 1.56,
			'MgMg': 2.0,
			'MgO': 1.83
		}

		self._atomic_masses = { #masses of the elements
			'C': 12,
			'H': 1,
			'O': 16,
			'N': 14,
			'Fe': 55.85,
			'Mg': 24.3
		}

		self._atom_colours = {
			'C': (34, 34, 34),
			'H': (255, 255, 255),
			'O': (255, 22, 0),
			'N': (22, 33, 255),
			'S': (225, 225, 48),
			'Ca': (61, 255, 0),
			'Fe': (221, 119, 0),
			'Mg': (0, 119, 0)
			}

		self._atom_radii = {
			'C': 0.67,
			'H': 0.53,
			'O': 0.48,
			'N': 0.56,
			'S': 1.00,
			'Ca': 1.80,
			'Fe': 1.40,
			'Mg': 1.50,
		}

		self.position = position	
		self.type = 'molecule'
		self.rotation = rotation
		self.scale = 400
		self.set_bonds()


	def __str__(self):
		'''
		Dunder method specifying what should be returned when the object is cast to a string.
		'''

		string = f'''
Molecule {self.name}:

Number of carbon atoms: {self.ncarbons}
Number of hydrogen atoms: {self._nelement("H")}

Coordinates (angstrom):
'''
		coord_string = [f'\t{a: ^4} {c[0]: >10.5f} {c[1]: >10.5f} {c[2]: >10.5f}\n' for a, c in zip(self.atom_types, self.coords)]
		string += ''.join(coord_string)
		string += f'\n\tCOM  {self.center_of_mass[0]: >10.5f} {self.center_of_mass[1]: >10.5f} {self.center_of_mass[2]: >10.5f}\t#Center of mass'

		return string


	def center(self):
		self.coords -= self.center_of_mass

	def _load_xyz(self, file):
		'''
		Method used to load xyz files

		file - string representing the path to the xyz file

		returns nothing, modifies self
		'''

		self.name = file.split('/')[-1].strip('.xyz')
		self.coords = np.loadtxt(file, skiprows=2, usecols=(1,2,3), dtype=float)
		self.atom_types = np.loadtxt(file, skiprows=2, usecols=0, dtype=str)


	@property 
	def ncarbons(self):
		'''	
		Property method that returns the number of carbon atoms in the molecule
		'''

		#use self._nelement method to return nmb of C
		return self._nelement('C')


	def _nelement(self, element):
		'''
		Method that returns the count of element in the molecule.
		We separate this from self._load_xyz in case we add support for other 
		file types later on (such as .sdf or .pdb).

		element - string of element to be counted

		returns integer count of element in molecule
		'''

		#cast self.atom_types to a list and count the number of 'C'
		return self.atom_types.tolist().count(element)


	def bond_angle(self, a1, a2, a3, in_degrees=False):
		'''
		Method that returns the a1--a2--a3 bond angle

		a1 - integer index of atom 1
		a2 - integer index of atom 2
		a3 - integer index of atom 3
		in_degrees - boolean specifying whether to return angle in degrees or radians
					 set to True for degrees or False for radians

		returns float a1--a2--a3 bond angle
		'''

		#get the two bond vectors:
		#visual example: 
		'''	
		a1	
		^	
		｜
		｜ u	
		｜ 
		｜	  v
		a2－－－－－> a3

		'''

		u = self.coords[a1] - self.coords[a2]
		v = self.coords[a3] - self.coords[a2]

		#We know that cos(theta) = u @ v / (|u| * |v|)
		#function to calculate the magnitude of a vector
		mag = lambda x: np.sqrt(x @ x)

		#return the angle. If in_degrees is True multiply it by 180/pi, else multiply by 1
		return np.arccos((u @ v) / (mag(u) * mag(v))) * (1, 180/math.pi)[in_degrees]


	def distance(self, a1, a2):
		'''
		Method that returns the euclidean distance between two atoms

		a1 - integer index of atom 1 or coordinates
		a2 - integer index of atom 2 or coordinates

		returns float distance between a1 and a2
		'''
		if type(a1) is int:
			a1 = self.coords[a1]
		if type(a2) is int:
			a2 = self.coords[a2]

		return euclidean(a1, a2)


	def isbonded(self, a1, a2):
		'''
		Method that returns whether two atoms are considered to be bonded.

		a1 - integer index of atom 1
		a2 - integer index of atom 2

		returns boolean 
		'''

		#get the atom_types of a1 and a2
		#sort the elements so that we can use self._bond_len_thresh which 
		#expects alphabetically sorted keys
		e1, e2 = sorted((self.atom_types[a1], self.atom_types[a2]))

		#get the distance between a1 and a2 and compare them to the e1-e2 bond
		#length threshold set in self.__init__

		return self.distance(a1, a2) < self._bond_len_thresh[e1 + e2]
		# return self.distance(a1, a2) < (rad[elem[a1]] + rad[elem[a2]]) * 1.1


	@property 
	def masses(self):
		'''
		Property method that returns the masses of the atoms in the molecule
		as an np.array
		'''

		return np.asarray([self._atomic_masses[e] for e in self.atom_types])


	@property
	def center_of_mass(self):
		'''
		Property method that returns the center of mass:
		COM = 1/M * sum[i over atoms](mi * ri) with M sum of all masses, mi
		mass of atom i, ri coords of atom i

		returns array of coordinates of center of mass
		'''

		#we simply multiply the (n, ) row vector self.masses by the (n,3) 
		#self.coords matrix. Matrix multiplication and dividing by M gives 
		#us the (1,3) coordinates of the COM.
		return self.masses @ self.coords / np.sum(self.masses) 


	def set_bonds(self):
		self.bonds = [self.get_bonded_atoms(a) for a in range(len(self.atom_types))]


	def get_bonded_atoms(self, a, element='any'):
		'''
		Method that returns the indices of atoms bonded to atom a

		a - integer index of atom a

		returns list of integers as indices of atoms bonded to a
		'''

		#Let i loop over the atoms and return i if it is bonded to a and i is not a 
		if element == 'any':
			return [i for i in range(len(self.atom_types)) if self.isbonded(a, i) and i != a]
		return [i for i in range(len(self.atom_types)) if self.isbonded(a, i) and i != a and self.atom_types[i] == element]


	def nbonds(self, a, element='any'):
		'''
		Method that returns the number of bonds to atom a.

		a - integer index of atom a

		returns integer number of bonds to a
		'''

		#First make a list of booleans (in python equivalent to 0 and 1) based on
		#self.isbonded between atom a and the rest of the atoms.
		#Then return the sum of the list.
		if element == 'any':
			return sum([self.isbonded(a, i) for i in range(len(self.atom_types)) if i != a])
		return sum([self.isbonded(a, i) for i in range(len(self.atom_types)) if i != a and self.atom_types[i] == element])


	def add_hydrogens(self, bond_len=1.1):
		'''
		Method that adds missing hydrogens to the molecule based on CCH bond angles
		of 120 deg and bond length for C-H of 1.1 Angstrom.
		Here we shall assume that the molecule is an aromatic hydrocarbon with CCH bond angles of around 120 deg.
		'''

		#Loop over the atoms
		new_H_coords = []
		for i, a in enumerate(self.atom_types):
			#check if atom is carbon and check if there are only two bonds to
			#the atom. Because the molecule if a PAH we know that carbons with only two
			#C-C bonds must also get a hydrogen.
			
			if a == 'C' and self.nbonds(i, element='C') == 2:
				#get the atoms bonded to atom a
				'''
				The strategy we will use to find the new position of the hydrogen is to first get the bond
				vectors of the C-C-C frame originating from C2.
				We then add the two vectors together to get the additive vector u + v
				Normalizing this vector and multiplying by the desired bond length for the C-H bond
				gives us the C-H bond vector in the opposite direction (because we intersected the inner
				C-C-C bond angle (~120deg), instead of the outer one (~240deg)). Therefore we subtract this
				vector from the coordinates of C2 to get the desired coordinates for the H-atom.

						   ^
						  /
						 /
						/
				C1	   /
				^	  /
				｜   /
			   u｜  /v+u	
				｜ /
				｜╱  
				C2－－－－－> C3
					  v

				Using this method we are also able to add hydrogens to non-planar aromatic hydrocarbons such as 
				hexabenzocoronene (molecules/hexabenzocoronene.xyz) or PAH in any orientation (need not be 
				aligned to the z-plane).
				'''

				#get the neighbours of the current carbon atom
				bonds = self.get_bonded_atoms(i, element='C')
				#calculate u an v
				v = self.coords[bonds[0]] - self.coords[i] 
				u = self.coords[bonds[1]] - self.coords[i] 
				#calculate and normalize (v + u)
				k = (v + u)/(np.sqrt((v + u) @ (v + u)))
				#multiply the unit vector by the bond length
				k *= bond_len
				# k = (v + u)/((self.distance(bonds[0], i) + self.distance(bonds[1], i))/2) * bond_len

				#append the list of new hydrogens
				new_H_coords.append(self.coords[i] - np.expand_dims(k, 0))
		#add all hydrogen atoms to the molecule. We do this last, so that new hydrogens do not interfere with
		#the above loop.
		[self.add_atom('H', c) for c in new_H_coords]


	def remove_hydrogens(self):
		'''
		Method that removes hydrogens from the molecule.
		'''
		#we loop over all atoms in reverse, so that we do not get complications
		#with out of bounds errors
		for i, a in zip(range(len(self.atom_types))[::-1], self.atom_types[::-1]):
			if a == 'H':
				self.atom_types = np.delete(self.atom_types, i)
				self.coords = np.delete(self.coords, i, axis=0)

		self.set_bonds()




	def add_atom(self, element, coords):
		'''
		Method that adds an atom to the molecule

		element - string element type of the atom
		coords - np.array with shape (1,3) representing the coordinates of the atom
		'''

		self.atom_types = np.append(self.atom_types, element)
		self.coords = np.append(self.coords, coords, axis=0)
		self.set_bonds()


	def save_to_xyz(self, file):
		'''
		Method that saves the molecule to a file

		file - string path to new file
		'''

		with open(file, 'w+') as f:
			f.write(f'{len(self.atom_types)}\n')
			f.write('generated via python\n')
			for e, c in zip(self.atom_types, self.coords):
				f.write(f'{e} {c[0]:.5f} {c[1]:.5f} {c[2]:.5f}\n')


	def rotate(self, rotation):
		'''
		Method that rotates atom coordinates by rotation
		'''

		r = rotation[0]
		Rx = np.array(([	  1, 	  0,	   0],
					   [	  0, cos(r), -sin(r)],
					   [      0, sin(r),  cos(r)]))

		r = rotation[1]
		Ry = np.array(([ cos(r),  	   0, sin(r)],
					   [ 	  0, 	   1,	   0],
					   [-sin(r), 	   0, cos(r)]))

		r = rotation[2]
		Rz = np.array(([ cos(r), -sin(r), 	   0],
					   [ sin(r),  cos(r), 	   0],
					   [ 	  0, 	   0, 	   1]))

		self.coords = np.asarray([Rx @ Ry @ Rz @ p for p in (self.coords - self.center_of_mass)])
