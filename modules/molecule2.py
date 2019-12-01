import numpy as np
from scipy.spatial.distance import euclidean
from math import cos, sin, pi
import os


class Atom:
	def __init__(self, element, coords):
		self.coords = coords
		self.element = element
		self.bonds = []
		self.bond_orders = {}

		self.mass = { #masses of the elements
			'C': 12,
			'H': 1,
			'O': 16,
			'N': 14,
			'Fe': 55.85,
			'Mg': 24.3,
			'P': 31,
			'Cl': 35.5,
			'S': 32.02,
			'Na': 23,
		}[element]

		self.colour = {
			'C': (34, 34, 34),
			'H': (255, 255, 255),
			'O': (255, 22, 0),
			'N': (22, 33, 255),
			'S': (225, 225, 48),
			'Ca': (61, 255, 0),
			'Fe': (221, 119, 0),
			'Mg': (0, 119, 0),
			'P': (255, 153, 0),
			'Cl': (31, 240, 31),
			'Na': (119, 0, 255),
			}[element]

		self.radius = {
			'Ac': 1.88,
			'Ag': 1.59,
			'Al': 1.35,
			'Am': 1.51,
			'As': 1.21,
			'Au': 1.50,
			'B': 0.83,
			'Ba': 1.34,
			'Be': 0.35,
			'Bi': 1.54,
			'Br': 0.68,
			'C': 0.68,
			'Ca': 0.99,
			'Cd': 1.69,
			'Ce': 1.83,
			'Cl': 0.99,
			'Co': 1.33,
			'Cr': 1.35,
			'Cs': 1.67,
			'Cu': 1.52,
			'D': 0.23,
			'Dy': 1.75,
			'Er': 1.73,
			'Eu': 1.99,
			'F': 0.64,
			'Fe': 1.34,
			'Ga': 1.22,
			'Gd': 1.79,
			'Ge': 1.17,
			'H': 0.23,
			'Hf': 1.57,
			'O': 0.68,
			'N': 0.68,
			'S': 1.02,

			
			'Mg': 1.50,
			'P': 1.10,
			
			'Na': 1.80,
		}[element]

		self.max_valence = {
			'C': 4,
			'H': 1,
			'O': 2,
			'N': 3,
			'Mg': 2,
			'P': 6,
			'Cl': 1,
			'S': 2,
			'Na': 1,
			'Fe': 2,
		}[element]


	def distance_to(self, p):
		'''
		Method that returns the distance to a point p. If p is of class Atom, use p's coords.

		p - float or Atom type.

		return float
		'''

		if type(p) is Atom:
			p = p.coords

		return euclidean(self.coords, p)


	def fix_overbonding(self):
		'''
		Method that attempts to fix overbonding.
		It does so by elminating the longest bond.
		'''

		while len(self.bonds) > self.max_valence:
			#sort bonds by distance
			self.bonds = sorted(self.bonds, key=lambda b: self.distance_to(b))
			#remove last item from sorted list
			b = self.bonds[-1]
			del(self.bonds[-1])
			#remove self from the neighbours bonds list

			b.bonds.remove(self)


	def set_bonds(self, molecule):
		'''
		Method that determines bonds based on  
		'''

		for atom in molecule.atoms:
			if not atom == self:
				if self.distance_to(atom) < self.radius + atom.radius + 0.4:
					self.bonds.append(atom)

		self.fix_overbonding()

		for b in self.bonds:
			self.bond_orders[b] = 1


	def get_bonds_by_elements(self, elements, blacklist=False):
		'''
		Method that returns the bonds to this element based
		on the elements of the other atoms.

		elements - list of strings specifying the elements to check
		blacklist - set to True to exclude elements, or to False to 
					include the elements

		returns list of Atom objects
		'''

		if blacklist:
			return [a for a in self.bonds if a.element not in elements]
		return [a for a in self.bonds if a.element in elements]


	def __repr__(self):
		return f'Atom({self.element}, {self.coords})'


	def is_unsaturated(self):
		return sum(self.bond_orders.values()) < self.max_valence

	def is_saturated(self):
		return sum(self.bond_orders.values()) == self.max_valence




class Molecule:
	'''
	Class representation of a molecule
	'''

	def __init__(self, position=[0.,0.,0.], rotation=[0.,0.,0.], molecule_file=None):
		self._warning_level = 1

		self.position = position	
		self.type = 'molecule'
		self.rotation = rotation
		self.scale = 400

		self.max_valence = {
			'C': 4,
			'H': 1,
			'O': 2,
			'N': 3,
			'Mg': 2,
			'P': 6,
			'Cl': 1,
			'S': 2,
			'Na': 1,
		}

		#check if file ends with xyz and try to load it
		if molecule_file is not None and molecule_file.endswith('.xyz'):
			self._load_xyz(molecule_file)

		if molecule_file is not None and molecule_file.endswith('.pcp'):
			self._load_from_pubchem(molecule_file.strip('.pcp'))


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
		for a in self.atoms:
			a.coords -= self.center_of_mass


	def _load_xyz(self, file):
		'''
		Method used to load xyz files

		file - string representing the path to the xyz file

		returns nothing, modifies self
		'''

		self.name = file.split('/')[-1].strip('.xyz')
		coords = np.loadtxt(file, skiprows=2, usecols=(1,2,3), dtype=float)
		elements = np.loadtxt(file, skiprows=2, usecols=0, dtype=str)

		self.atoms = []
		[self.atoms.append(Atom(elements[i], coords[i])) for i in range(len(coords))]

		self.set_bonds()


	def _load_from_pubchem(self, name):
		'''
		Method used to load data from pubchem

		name - name of compound to search

		returns nothing
		'''

		import pubchempy as pcp

		try:
			name = int(name)
		except:
			pass

		record_type = '3d'
		mol = pcp.get_compounds(name, ('name', 'cid')[type(name) is int], record_type=record_type)

		if len(mol) == 0:
			print(f'Could not find 3d structure of {name}... Attempting to find 2d structure...')
			record_type = '2d'
			mol = pcp.get_compounds(name, ('name', 'cid')[type(name) is int], record_type=record_type)

		if len(mol) == 0:
			print(f'No structural data found for {name}')

		else:
			mol = mol[0]

			coords = np.asarray([[a.x, a.y, a.z] for a in mol.atoms])
			coords = np.where(coords == None, 0, coords).astype(float)
			elements = np.asarray([a.element for a in mol.atoms])
			self.name = name

			self.atoms = []
			[self.atoms.append(Atom(elements[i], coords[i])) for i in range(len(coords))]
			if record_type == '3d':
				self.save_to_xyz(os.getcwd() + rf'\Molecules\{name.lower()}.xyz')

			self.set_bonds()


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
		return [a.element for a in self.atoms].count(element)


	def bond_angle(self, a1, a2, a3, in_degrees=False):
		'''
		Method that returns the a1--a2--a3 bond angle

		a1 - atom 1
		a2 - atom 2
		a3 - atom 3
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

		u = a1.coords - a2.coords
		v = a3.coords - a2.coords

		#We know that cos(theta) = u @ v / (|u| * |v|)
		#function to calculate the magnitude of a vector
		mag = lambda x: np.sqrt(x @ x)

		#return the angle. If in_degrees is True multiply it by 180/pi, else multiply by 1
		return np.arccos((u @ v) / (mag(u) * mag(v))) * (1, 180/pi)[in_degrees]


	def set_HA_valence(self):
		for a in self.atoms:
			a.HA_valence = len(a.get_bonds_by_elements(['H'], blacklist=True))


	def set_hybridisation(self):
		for a in self.atoms:
			if a.element == 'C':
				if a.HA_valence == 4:
					a.hybridisation = 'sp3'

				if a.HA_valence == 3:
					bonds = a.get_bonds_by_elements(['H'], blacklist=True)
					bond_angles = self.bond_angle(bonds[0], a, bonds[1], in_degrees=True)
					bond_angles += self.bond_angle(bonds[0], a, bonds[2], in_degrees=True)
					bond_angles += self.bond_angle(bonds[1], a, bonds[2], in_degrees=True)
					bond_angle = bond_angles/3
					if bond_angle > 115.0:
						a.hybridisation = 'sp2'
					else:
						a.hybridisation = 'sp3'

				if a.HA_valence == 2:
					bonds = a.get_bonds_by_elements(['H'], blacklist=True)
					bond_angle = self.bond_angle(bonds[0], a, bonds[1], in_degrees=True)
					if bond_angle > 160.0:
						a.hybridisation = 'sp'
					elif bond_angle > 115.0:
						a.hybridisation = 'sp2'
					else:
						a.hybridisation = 'sp3'

				if a.HA_valence == 1: 
					a.hybridisation = 'sp3'
			else: a.hybridisation = 'sp3'


	def distance(self, a1, a2):
		'''
		Method that returns the euclidean distance between two atoms

		a1 - atom 1 or coordinates
		a2 - atom 2 or coordinates

		returns float distance between a1 and a2
		'''
		if type(a1) is Atom:
			a1 = a1.coords
		if type(a2) is Atom:
			a2 = a2.coords

		return euclidean(a1, a2)


	def isbonded(self, a1, a2):
		'''
		Method that returns whether two atoms are considered to be bonded.

		a1 - atom 1
		a2 - atom 2

		returns boolean 
		'''

		return a2 in a1.bonds


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
		masses = np.asarray([a.mass for a in self.atoms])
		return masses @ np.asarray([a.coords for a in self.atoms]) / np.sum(masses) 


	def set_bonds(self):
		for a in self.atoms:
			a.bonds = []
			a.bond_orders = {}

		self.guess_bond_order_iters = 0
		[a.set_bonds(self) for a in self.atoms]

		self.fix_overbonding()
		self.set_HA_valence()
		self.set_hybridisation()
		self.guess_bond_orders()


	def reset_bond_orders(self):
		for a in self.atoms:
			for key in a.bond_orders.keys():
				a.bond_orders[key] = 1

	def guess_bond_orders(self):
		'''
		Method that guesses the bond orders of the molecule.
		
		Current strategy:
		- Sort elements from low valence to high valence (H < O < N < C, etc..) 
		  and loops over the elements.
			- Collect every atom of the element and checks its bond saturation.
			- If the atom is not saturated, loop over the atoms it is bonded to.
				- Check the saturation of the bonded atom. If the bonded atom 
				  is also not saturated, increase the bond order to that bond.
				- Terminate the loop if the current atom is saturated.

		Problems:
		- Sometimes all but one atom is saturated. For example with chlorophyll.

		Possible Solution:
		- Check at the end if there are unsaturated atoms. Redo the calculation but in different orders
		'''

		self.reset_bond_orders()

		#sort the elements by valence
		valences = list(self.max_valence.items())
		valences = sorted(valences, key=lambda x: x[1])

		for el, _ in valences:
			#get all atoms of element el
			atoms = self.get_by_element(el).copy()
			np.random.shuffle(atoms)
			#loop over the atoms
			for a in atoms:
				#check if a is saturated
				if a.is_unsaturated():
					#if not, get its neighbours
					neighbours = np.copy(self.get_bonded_atoms(a))
					neighbours = sorted(neighbours, key=lambda x: a.distance_to(x))
					neighbour_hybrids = [x.hybridisation for x in neighbours]
					#loop over the neighbours

					# if a.hybridisation == 'sp' and 'sp' in neighbour_hybrids:
					# 	x = neighbour_hybrids.index('sp')
					# 	a.bond_orders[neighbours[x]] += 2
					# 	neighbour.bond_orders[a] += 2

					if a.is_unsaturated():
						for neighbour in neighbours:
							#check if the neighbour is also unsaturated and whether a has 
							#become saturated  in the meantime
							if neighbour.is_unsaturated() and a.is_unsaturated():
								a.bond_orders[neighbour] += 1
								neighbour.bond_orders[a] += 1


		#give warnings if necessary
		mbo = sum([a.is_saturated() for a in self.get_by_element('C')]) - len(self.get_by_element('C'))
		if self._warning_level >= 1:
			if mbo < 0:
				print(f'Molecule.guess_bond_orders: Bond order guessing was not succesful. Unsaturated atoms: {abs(mbo)} (iteration {self.guess_bond_order_iters})')
			else:
				print('Molecule.guess_bond_orders: Bond orders seem fine')
		
		if mbo < 0 and self.guess_bond_order_iters < self.natoms:
			self.guess_bond_order_iters += 1
			self.guess_bond_orders()
			


	@property
	def natoms(self):
		return len(self.atoms)


	def is_saturated(self, a):
		return self.get_saturation(a) == self._atom_valence[self.atom_types[a]]


	def is_unsaturated(self, a):
		return self.get_saturation(a) < self._atom_valence[self.atom_types[a]]


	def fix_overbonding(self):
		[a.fix_overbonding() for a in self.atoms]


	def get_saturation(self, a):
		'''
		Method that returns the valence saturation of an atom.

		a - integer index of atom

		return integer number of bonds
		'''

		return sum(self.bond_orders[a])


	def get_by_element(self, element):
		'''
		Method that returns a list of atoms corresponding to the atoms in the 
		molecule of a certain element.

		element - string specifying the element

		returns list of indices
		'''

		return [a for a in self.atoms if a.element == element]


	def get_bonded_atoms(self, a, element='any'):
		'''
		Method that returns the indices of atoms bonded to atom a

		a - integer index of atom a

		returns list of integers as indices of atoms bonded to a
		'''

		if element == 'any':
			return a.bonds
		else:
			return [b for b in a.bonds if b.element == element]


	def nbonds(self, a, element='any'):
		'''
		Method that returns the number of bonds to atom a.

		a - integer index of atom a

		returns integer number of bonds to a
		'''

		return len(get_bonded_atoms(a, element))


	def add_hydrogens(self, bond_len=1.1):
		'''
		Method that adds missing hydrogens to the molecule based on CCH bond angles
		of 120 deg and bond length for C-H of 1.1 Angstrom.
		Here we shall assume that the molecule is an aromatic hydrocarbon with CCH bond angles of around 120 deg.
		'''

		#Loop over the atoms
		new_H_coords = []
		for i, a in enumerate(self.atoms):
			#check if atom is carbon and check if there are only two bonds to
			#the atom. Because the molecule if a PAH we know that carbons with only two
			#C-C bonds must also get a hydrogen.
			
			if a.element == 'C' and len(a.get_bonds_by_elements('C')) == 2:
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
				bonds = a.get_bonds_by_elements('C')
				#calculate u an v
				v = bonds[0].coords - a.coords
				u = bonds[1].coords - a.coords
				#calculate and normalize (v + u)
				k = (v + u)/(np.sqrt((v + u) @ (v + u)))
				#multiply the unit vector by the bond length
				k *= bond_len
				# k = (v + u)/((self.distance(bonds[0], i) + self.distance(bonds[1], i))/2) * bond_len

				#append the list of new hydrogens
				new_H_coords.append(a.coords - np.expand_dims(k, 0))

		#add all hydrogen atoms to the molecule. We do this last, so that new hydrogens do not interfere with
		#the above loop.
		[self.add_atom('H', c, set_bonds=False) for c in new_H_coords]

		self.set_bonds()

	def remove_by_element(self, element):
		'''
		Method that removes hydrogens from the molecule.

		element - string element type of the atom
		'''

		[self.atoms.remove(a) for a in self.get_by_element(element)]
		self.set_bonds()


	def remove_atom(self, a):
		'''
		Method that removes an atom from the molecule.

		a - integer index of atom a
		'''

		self.atom_types = np.delete(self.atom_types, a)
		self.coords = np.delete(self.coords, a, axis=0)


	def add_atom(self, element, coords, set_bonds=True):
		'''
		Method that adds an atom to the molecule

		element - string element type of the atom
		coords - np.array with shape (1,3) representing the coordinates of the atom
		'''

		self.atoms.append(Atom(element, coords.flatten()))


	def save_to_xyz(self, file):
		'''
		Method that saves the molecule to a file

		file - string path to new file
		'''

		with open(file, 'w+') as f:
			f.write(f'{self.natoms}\n')
			f.write('generated via python\n')
			for a in self.atoms:
				f.write(f'{a.element} {a.coords[0]:.5f} {a.coords[1]:.5f} {a.coords[2]:.5f}\n')


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

		for a in self.atoms:
			a.coords = Rx @ Ry @ Rz @ a.coords

