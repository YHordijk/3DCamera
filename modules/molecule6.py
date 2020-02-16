import numpy as np
from scipy.spatial.distance import euclidean, sqeuclidean
from math import cos, sin, pi, atan2, acos, exp
import os, json, time
import modules.basisset6 as bs
import modules.utils as utils
import periodictable as pt
import networkx as nx




class Atom:
	def __init__(self, element, coords):
		self.coords = coords
		self.bonds = []
		self.bond_orders = {}
		self.selected = False

		try:
			el = pt.elements[element]
		except:
			try:
				el = pt.elements.symbol(element)
			except:
				try:
					el = pt.elements.name(element)
				except:
					utils.message('Error: Could not parse element {element}.', colour='red')

		self.ionisation_energy = np.genfromtxt('modules\\ionisation_energies', usecols=(1,2), missing_values='', delimiter=';')[el.number]

		self.symbol = el.symbol
		self.name = el.name
		self.atomic_number = el.number
		self.mass = el.mass
		self.radius = el.covalent_radius
		self.electro_negativity = {
			1:2.2,
			2:0,
			3:0.98,
			4:1.57,
			5:2.04,
			6:2.55,
			7:3.04,
			8:3.44,
			9:3.98,
			10:0,
			11:0.93,
			12:1.31,
			13:1.61,
			14:1.9,
			15:2.19,
			16:2.58,
			17:3.16,
			18:0,
			19:0.82,
			20:1,
			21:1.36,
			22:1.54,
			23:1.63,
			24:1.66,
			25:1.55,
			26:1.83,
			27:1.88,
			28:1.91,
			29:1.9,
			30:1.65,
			31:1.81,
			32:2.01,
			33:2.18,
			34:2.55,
			35:2.96,
			36:3,
			37:0.82,
			38:0.95,
			39:1.22,
			40:1.33,
			41:1.6,
			42:2.16,
			43:1.9,
			44:2.2,
			45:2.28,
			46:2.2,
			47:1.93,
			48:1.69,
			49:1.78,
			50:1.96,
			51:2.05,
			52:2.1,
			53:2.66,
			54:2.6,
			55:0.79,
			56:0.89,
			57:1.1,
			58:1.12,
			59:1.13,
			60:1.14,
			61:1.13,
			62:1.17,
			63:1.2,
			64:1.2,
			65:1.22,
			66:1.23,
			67:1.24,
			68:1.24,
			69:1.25,
			70:1.1,
			71:1.27,
			72:1.3,
			73:1.5,
			74:2.36,
			75:1.9,
			76:2.2,
			77:2.2,
			78:2.28,
			79:2.54,
			80:2,
			81:1.62,
			82:2.33,
			83:2.02,
			84:2,
			85:2.2,
			86:0,
			87:0.7,
			88:0.89,
			89:1.1,
			90:1.3,
			91:1.5,
			92:1.38,
			93:1.36,
			94:1.28,
			95:1.3,
			96:1.3,
			97:1.3,
			98:1.3,
			99:1.3,
			100:1.3,
			101:1.3,
			102:1.3,
			103:00,
			104:0,
			105:0,
			106:0,
			107:0,
			108:0,
			109:0,
			110:0,
			111:0,
			112:0,
			113:0,
			114:0,
			115:0,
			116:0,
			117:0,
			118:0,
		}[self.atomic_number]

		try:
			self.GMP_electro_negativity = {
				1: 0.89,
				3: 0.97,
				4: 1.47,
				5: 1.6,
				6: 2,
				7: 2.61,
				8: 3.15,
				9: 3.98,
				11: 1.01,
				12: 1.23,
				13: 1.47,
				14: 1.58,
				15: 1.96,
				16: 2.35,
				17: 2.74,
				19: 0.91,
				20: 1.04,
				21: 1.2,
				22: 1.32,
				23: 1.45,
				24: 1.56,
				25: 1.6,
				26: 1.64,
				27: 1.7,
				28: 1.75,
				29: 1.75,
				30: 1.66,
				31: 1.82,
				32: 1.51,
				33: 2.23,
				34: 2.51,
				35: 2.58,
				37: 0.89,
				38: 0.99,
				39: 1.11,
				40: 1.22,
				41: 1.23,
				42: 1.3,
				44: 1.42,
				45: 1.54,
				46: 1.35,
				47: 1.42,
				48: 1.46,
				49: 1.49,
				50: 1.72,
				51: 1.72,
				52: 2.72,
				53: 2.38,
				55: 0.86,
				56: 0.97,
				57: 1.08,
				58: 1.08,
				59: 1.07,
				60: 1.07,
				62: 1.07,
				63: 1.01,
				64: 1.11,
				65: 1.1,
				66: 1.1,
				67: 1.1,
				68: 1.11,
				69: 1.11,
				70: 1.06,
				71: 1.14,
				72: 1.23,
				73: 1.33,
				74: 1.4,
				75: 1.46,
				77: 1.55,
				80: 1.44,
				81: 1.44,
				82: 1.55,
				83: 1.67,
				90: 1.11,
				92: 1.22,
			}[self.atomic_number]
		except:
			self.GMP_electro_negativity = 0


		try:
			self.colour = self.draw_colour = {
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
				}[self.symbol]
		except:
			utils.message('Error: No default colour found for {self.symbol}. Defaulting to (0,0,0).', colour='red')
			self.colour = (0,0,0)

		try:
			self.max_valence = {
				'C': 4,
				'H': 1,
				'O': 2,
				'N': 3,
				'Mg': 2,
				'P': 6,
				'Cl': 1,
				'S': 6,
				'Na': 1,
				'Fe': 2,
				}[self.symbol]
		except:
			utils.message('Error: No max_valence found for {self.symbol}. Defaulting to 1.', colour='red')
			self.max_valence = 1


	def bond_dist_to(self, a):
		'''
		Method that returns the number of bonds to atom a.
		'''

		bdist = 0
		atoms = [self]
		while a not in atoms:
			bdist += 1
			new_atoms = []
			for atom in atoms:
				new_atoms += atom.bonds
			atoms = new_atoms
			if bdist > 10:
				break

		return bdist

	def distance_to(self, p):
		'''
		Method that returns the distance to a point p. If p is of class Atom, use p's coords.

		p - tuple of float or Atom type.

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
			try:
				del(b.bond_orders[self])
				del(self.bond_orders[b])
			except:
				pass

			del(self.bonds[-1])

			#remove self from the neighbours bonds list
			b.bonds.remove(self)


	@property 
	def valence(self):
		return sum([bo for bo in self.bond_orders.values()])


	@property 
	def penalty_score(self):
		return self.max_valence - self.valence


	def set_bonds(self, molecule):
		'''
		Method that determines bonds based on  
		'''

		for atom in molecule.atoms:
			if not atom == self:
				if self.distance_to(atom) < self.radius + atom.radius + 0.4:
					self.bonds.append(atom)

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
			return [a for a in self.bonds if a.symbol not in elements]
		return [a for a in self.bonds if a.symbol in elements]


	def __repr__(self):
		return f'Atom({self.symbol}, {self.coords})'


	def is_unsaturated(self):
		return sum(self.bond_orders.values()) < self.max_valence

	def is_saturated(self):
		return sum(self.bond_orders.values()) == self.max_valence


class Molecule:
	'''
	Class representation of a molecule
	'''

	def __init__(self, molecule_file=None, atoms=None, position=[0.,0.,0.], rotation=[0.,0.], 
					warning_level=1, scale=200, basis_set_type='STO-6G', repeat=1, repeat_vector=None):
		self._warning_level = warning_level

		self.position = position	
		self.type = 'molecule'
		self.rotation = rotation
		self.scale = scale
		self.basis_set_type = basis_set_type
		self.repeat = repeat
		self.repeat_vector = repeat_vector

		self._dens_pos = {}
		self._dens_colours = {}

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

		if atoms is not None:
			self.atoms = atoms
			self.name = 'Molecule'
			self._mol_load_finish()

		else:
			self._load_mol(molecule_file)


	def _load_mol(self, molecule_file):
		#get the file extension:
		if '.' in molecule_file:
			name, filetype = molecule_file.split('.')
		else:
			utils.message(f'Searching pubchem for {molecule_file}...')
			self._load_from_pubchem(molecule_file)
			return

		if filetype == 'pcp':
			utils.message(f'Searching pubchem for {name}...')
			self._load_from_pubchem(name)
			return

		else:
			#check if file exists:
			if not os.path.exists(molecule_file):
				utils.message(f'Error: File {molecule_file} does not exist. Searching pubchem ...', colour='red')
				self._load_from_pubchem(molecule_file)
				return
			else:
				if filetype == 'xyz':
					self._load_xyz(molecule_file)
					return
				if filetype == 'xyzb':
					self._load_xyzb(molecule_file)
					return


		# self._load_from_file()
		# #check if file ends with xyz and try to load it
		# if molecule_file is not None and molecule_file.endswith('.xyz'):
			

		# # elif molecule_file is not None and molecule_file.endswith('.pcp'):
		# # 	self._load_from_pubchem(molecule_file[0:-4])

		# else:
		# 	self._load_xyz(os.getcwd() + f'\\Molecules\\{molecule_file}.xyz')



	def __str__(self):
		'''
		Dunder method specifying what should be returned when the object is cast to a string.
		'''

		string = ''
		coord_string = [f'\t{a.symbol: ^4} {a.coords[0]: >10.5f} {a.coords[1]: >10.5f} {a.coords[2]: >10.5f}\n' for a in self.atoms]
		string += ''.join(coord_string)

		return string


	def electrostatic_potential(self, pos):
		from scipy.spatial.distance import cdist

		p = np.expand_dims(np.zeros(pos.shape[0]),1)
		for atom in self.atoms:
			a = np.expand_dims(atom.coords, 1).T
			r = cdist(pos, a) # calculate separation
			p -= atom.atomic_number / r

		return p


	def center(self, c=None):
		if c is None:
			c = self.center_of_mass.copy()
		if type(c) is Atom:
			c = c.coords.copy()

		for a in self.atoms:
			a.coords -= c

		self.position = (0,0,0)


	def reset_colours(self):
		for a in self.atoms:
			a.draw_colour = a.colour
			

	def _load_xyz(self, file):
		'''
		Method used to load xyz files

		file - string representing the path to the xyz file

		returns nothing, modifies self
		'''
		self.name = file.split('/')[-1].strip('.xyz').capitalize()
		self.name = self.name.split('\\')[-1].strip('.xyz').capitalize()

		coords = list(np.loadtxt(file, skiprows=2, usecols=(1,2,3), dtype=float))
		elements = list(np.loadtxt(file, skiprows=2, usecols=0, dtype=str))

		if 'VEC1' in elements:
			i = elements.index('VEC1')
			self.repeat_vector = coords[i]
			del(elements[i])
			del(coords[i])

		self.atoms = []
		for r in range(self.repeat):
			for e, c in zip(elements, coords):
				if r == 0:
					self.atoms.append(Atom(e, c))
				else:
					if not type(self.repeat_vector) is np.ndarray:
						utils.message('Error: Please supply repeat vector.', colour='red')
					else:
						self.atoms.append(Atom(e, c+(r)*self.repeat_vector))

		self._mol_load_finish()


	def _load_xyzb(self, file):
		'''
		Method used to load xyz files

		file - string representing the path to the xyz file

		returns nothing, modifies self
		'''
		self.name = file.split('/')[-1].strip('.xyzb').capitalize()
		self.name = self.name.split('\\')[-1].strip('.xyzb').capitalize()

		coords = list(np.loadtxt(file, skiprows=2, usecols=(1,2,3), dtype=float))
		elements = list(np.loadtxt(file, skiprows=2, usecols=0, dtype=str))

		if 'VEC1' in elements:
			i = elements.index('VEC1')
			self.repeat_vector = coords[i]
			del(elements[i])
			del(coords[i])

		self.atoms = []
		for r in range(self.repeat):
			for e, c in zip(elements, coords):
				if r == 0:
					self.atoms.append(Atom(e, c))
				else:
					if not type(self.repeat_vector) is np.ndarray:
						utils.message('Error: Please supply repeat vector.', colour='red')
					else:
						self.atoms.append(Atom(e, c+(r)*self.repeat_vector))

		self._mol_load_finish()


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
			utils.message('Error: Could not find 3d structure of {name}... Attempting to find 2d structure...', colour='red')
			record_type = '2d'
			mol = pcp.get_compounds(name, ('name', 'cid')[type(name) is int], record_type=record_type)

		if len(mol) == 0:
			utils.message('Error: No structural data found for {name}.', colour='red')

		else:
			mol = mol[0]

			coords = np.asarray([[a.x, a.y, a.z] for a in mol.atoms])
			coords = np.where(coords == None, 0, coords).astype(float)
			elements = np.asarray([a.element for a in mol.atoms])
			self.name = name.capitalize()

			self.atoms = []
			[self.atoms.append(Atom(elements[i], coords[i])) for i in range(len(coords))]
			if record_type == '3d':
				self.save(os.getcwd() + rf'\Molecules\{name.lower()}.xyz')

			self._mol_load_finish()


	def _mol_load_finish(self):
		'''
		Method that is called by both xyz and pubchem loading of molecules
		'''

		self.center()

		self.set_bonds()
		self.detect_rings()

		self._load_basis_set()

		utils.message(f'Succesfully loaded {self.name}.', colour='green')



	def _load_basis_set(self):
		self.basis = bs.Basis(self, self.basis_set_type)


	def get_orb_density(self, p):
		'''
		Method that returns the density value at a point based on the selected basis set
		'''

		return self.basis_set(p)


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
		return [a.symbol for a in self.atoms].count(element)


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


	def torsion_angle(self, a1, a2, a3, a4, in_degrees=False):
		'''
		Method that returns the torsion angle or dihedral angle of the 
		a1 -- a2 -- a3 and a2 -- a3 -- a4 planes.

		a - atom object
		in_degrees - boolean specifying whether to return angle in degrees or radians
					 set to True for degrees or False for radians

		returns float
		'''

		#using method provided on https://en.wikipedia.org/wiki/Dihedral_angle
		b1 = a2.coords - a1.coords
		b2 = a3.coords - a2.coords
		b3 = a4.coords - a3.coords

		return atan2(np.dot(np.cross(np.cross(b1, b2), np.cross(b2, b3)), b2/np.linalg.norm(b2)), np.dot(np.cross(b1, b2), np.cross(b2, b3)))


	def set_HA_valence(self):
		for a in self.atoms:
			a.HA_valence = len(a.get_bonds_by_elements(['H'], blacklist=True))


	def set_hybridisation(self):
		for a in self.atoms:
			c = len(self.get_bonded_atoms(a))
			if a.max_valence == 1:
				a.hybridisation = 0

			elif a.max_valence == 2:
				if c == 2:
					a.hybridisation = 3
				elif c == 1:
					a. hybridisation = 2

			elif a.max_valence == 3:
				if c == 3:
					a.hybridisation = 3
				elif c == 2:
					a.hybridisation = 2
				elif c == 1:
					a.hybridisation = 1

			elif a.max_valence == 4:
				if c == 4:
					a.hybridisation = 3
				elif c == 3:
					a.hybridisation = 2
				elif c == 2:
					a.hybridisation = 1

		for a in self.atoms:
			if a.symbol == 'O' or a.symbol == 'N' and a.hybridisation == 3:
				if any([a2.hybridisation == 2 for a2 in a.bonds]):
					a.hybridisation = 2


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


	def detect_rings(self):
		'''
		Method that detects the rings in the molecule using
		networkx module. It also detects the type of ring
		(aliphatic (AL), aromatic (AR))
		'''

		#create a graph first
		g = nx.Graph()
		g.add_nodes_from(self.atoms)
		g.add_edges_from(list(self.get_unique_bonds()))
		#get the cycles
		cycles = nx.algorithms.cycles.minimum_cycle_basis(g)

		self.rings = []
		for cycle in cycles:
			#get some data on the atoms
			atom_sym = [a.symbol for a in cycle]
			carbon_hybrids = [a.hybridisation for a in cycle if a.symbol == 'C']
			nitrogen_bonds = [len(a.bonds) for a in cycle if a.symbol == 'N']

			#carbon contributes 1 electron, oxygen 2, nitrogen with 3 bonds 2, nitrogens with 2 bonds 1
			ne = carbon_hybrids.count(2) + atom_sym.count('O')*2 + nitrogen_bonds.count(3)*2 + nitrogen_bonds.count(2)

			#apply kekule rule and check if all carbons are sp2
			if all([h == 2 for h in carbon_hybrids]) and ne%4==2:
				self.rings.append((cycle, 'AR'))
			else:
				self.rings.append((cycle, 'AL'))

		#give the atoms ring properties
		for atom in self.atoms:
			atom.ring = None
			for cycle, typ in self.rings:
				if not atom.ring == 'AR':
					if atom in cycle:
						atom.ring = typ
					else:
						atom.ring = 'NO'


	def shake(self, intens=1):
		intens
		for a in self.atoms:
			a.coords += np.random.random(size=3) * intens *2 - intens
			

	def set_torsion_angle(self, a1, a2, a3, a4, theta):
		'''
		Method that sets the torsion angle between a1, a2, a3 and a4 to theta.
		'''
		self.rotate_bond(a2, a3, theta - self.torsion_angle(a1,a2,a3,a4))



	def rotate_bond(self, a1, a2, r):
		'''
		Method that rotates a bond between atoms a1 and a2 to r radians.

		a1, a2 - atom objects
		r - rotation in radians
		'''

		#get fragments
		frag1, frag2 = self.separate_on_bond(a1, a2)

		#get bond vector
		b = a1.coords - a2.coords

		Rx = np.array(([	  1, 	  0,	   0],
					   [	  0, cos(r), -sin(r)],
					   [      0, sin(r),  cos(r)]))

		self.align_bond_to_vector(a1,a2,(1,0,0))
		self.center(a1)

		for a in frag1:
			a.coords = Rx @ a.coords

		self.align_bond_to_vector(a1, a2, b)
		self.center()


	def align_bond_to_vector(self, a1, a2, b):
		'''
		Method that aligns a bond along a1 and a2 to the x-axis
		https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
		'''

		#get the bond vector
		a = a1.coords - a2.coords

		a = a/np.linalg.norm(a)
		b = b/np.linalg.norm(b)
		v = np.cross(a, b)
		c = np.dot(a, b)

		if not c == -1.0:
			skew = np.array([[0, -v[2], v[1]],
							 [v[2], 0, -v[0]],
							 [-v[1], v[0], 0]])

			R = np.eye(3) + skew + skew@skew * (1/(1+c))
			for a in self.atoms:
				a.coords = R @ a.coords


	def stretch_bond(self, a1, a2, r):
		'''
		Method that stretches a bond between atoms a1 and a2 to length r

		a1, a2 - atom objects
		r - desired bond length
		'''

		#get molecule fragments
		frag1, frag2 = self.separate_on_bond(a1, a2)

		#get the bond direction
		b = a1.coords - a2.coords
		b_norm = b/np.linalg.norm(b)

		#we only have to transpose one of the fragments
		for a in frag1:
			a.coords += - b + b_norm * r



	def separate_on_bond(self, a1, a2):
		'''
		Method that separates molecule on the specified bond. Will return
		the atoms that make up the two fragments. Fails if a1 and a2 are 
		part of the same ring or when a1 is not bonded to a2.

		a1, a2 - atom objects

		returns tuple of two lists
		'''

		#check if a1 and a2 are part of the same ring
		for ring, _ in self.rings:
			if a1 in ring and a2 in ring:
				return [], []

		#check if they are bonded:
		if not self.isbonded(a1, a2):
			return [], []

		#get atoms attached to a1:
		atoms1 = [a1]
		for a in atoms1:
			for b in a.bonds:
				if not b == a2:
					if not b in atoms1:
						atoms1.append(b)

		#get atoms attached to a2:
		atoms2 = [a2]
		for a in atoms2:
			for b in a.bonds:
				if not b == a1:
					if not b in atoms2:
						atoms2.append(b)

		return atoms1, atoms2


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


	def get_chiral_centres(self):
		for a in self.atoms:
			...


	def get_unique_atom_pairs(self):
		'''
		Generator method that yields all unique pairs of atoms in the molecule
		'''

		atoms = self.atoms
		for i in range(len(atoms)):
			a1 = atoms[i]
			for j in range(i+1, len(atoms)):
				a2 = atoms[j]
				yield sorted((a1,a2), key=lambda x: id(x))

	def get_unique_nonbonded_atom_pairs(self):
		'''
		Generator method that yields all unique pairs of atoms in the molecule
		'''

		atoms = self.atoms
		for i in range(len(atoms)):
			a1 = atoms[i]
			for j in range(i+1, len(atoms)):
				a2 = atoms[j]
				if not a2 in a1.bonds:
					yield sorted((a1,a2), key=lambda x: id(x))


	def get_unique_bonds(self):
		'''
		Generator method that yields all unique bonds in the molecule.
		'''

		prev_atoms = []

		for a1 in self.atoms:
			prev_atoms.append(a1)
			for a2 in a1.bonds:
				if not a2 in prev_atoms:
					yield (a1, a2)


	def get_unique_bond_angles(self, in_degrees=False):
		'''
		Generator method that yields all unique bond angles in the molecule along with the atoms over which the bond angle is calculated.
		'''

		prev_angles = []

		for a1 in self.atoms:
			for a2 in a1.bonds:
				for a3 in a2.bonds:
					sorted_atoms = sorted((a1,a2,a3), key=lambda x: id(x))
					if not sorted_atoms in prev_angles and len(set((a1, a2, a3))) == 3:
						prev_angles.append(sorted_atoms)
						yield (a1, a2, a3, self.bond_angle(a1, a2, a3, in_degrees=in_degrees))


	def get_unique_torsion_angles(self, in_degrees=False):
		'''
		Method that yields all unique torsion angles in the molecule along with the atoms over which the torsion angle is calculated.
		'''

		prev_angles = []
		for a1 in self.atoms:
			for a2 in a1.bonds:
				for a3 in a2.bonds:
					for a4 in a3.bonds:
						sorted_atoms = sorted((a1,a2,a3, a4), key=lambda x: id(x))
						if not sorted_atoms in prev_angles and len(set((a1, a2, a3, a4))) == 4:
							prev_angles.append(sorted_atoms)
							yield (a1, a2, a3, a4, self.torsion_angle(a1, a2, a3, a4, in_degrees))



	def get_rotatable_bonds(self):
		
		for a1, a2 in self.get_unique_bonds():
			if len(a1.bonds) > 1 and len(a2.bonds) > 1:
				yield (a1,a2)


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
		'''

		self.reset_bond_orders()

		#sort the elements by valence
		valences = list(self.max_valence.items())
		valences = sorted(valences, key=lambda x: x[1])

		for el, _ in valences:
			#get all atoms of element el
			atoms = self.get_by_element(el).copy()
			if self.guess_bond_order_iters > 0:
				np.random.shuffle(atoms)
			else:
				atoms = sorted(atoms, key=lambda x: x.hybridisation)
			#loop over the atoms
			for a in atoms:
				#check if a is saturated
				if a.is_unsaturated():
					#if not, get its neighbours
					neighbours = np.copy(self.get_bonded_atoms(a))
					# np.random.shuffle(neighbours)
					neighbours = sorted(neighbours, key=lambda x: a.distance_to(x))
					neighbour_hybrids = [x.hybridisation for x in neighbours]
					#loop over the neighbours

					if a.hybridisation == 'sp' and 'sp' in neighbour_hybrids:
						x = neighbour_hybrids.index('sp')
						a.bond_orders[neighbours[x]] = 3
						neighbour.bond_orders[a] = 3

					if a.is_unsaturated():
						for neighbour in neighbours:
							#check if the neighbour is also unsaturated and whether a has 
							#become saturated  in the meantime
							if neighbour.is_unsaturated() and a.is_unsaturated():
								a.bond_orders[neighbour] += 1
								neighbour.bond_orders[a] += 1


		#give warnings if necessary
		mbo = sum([a.is_saturated() for a in self.get_by_element('C')]) - len(self.get_by_element('C'))
		if self._warning_level == 2:
			if mbo < 0:
				utils.message(f'Error: Bond order guessing was not succesful. Unsaturated atoms: {abs(mbo)} (iteration {self.guess_bond_order_iters}).', colour='red')
			else:
				utils.message(f'Bond order guessing succesful after {self.guess_bond_order_iters+1} iterations', colour='green')

		elif self._warning_level == 1:
			if mbo == 0:
				utils.message(f'Bond order guessing succesful after {self.guess_bond_order_iters+1} iteration{"s"*(self.guess_bond_order_iters>0)}.', colour='green')
			elif self.guess_bond_order_iters == self.natoms:
				utils.message(f'Bond order guessing was not succesful. Unsaturated atoms: {abs(mbo)}.', colour='red')
		
		if mbo < 0 and self.guess_bond_order_iters < 5 * self.natoms:
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

	def get_elements(self):
		'''
		Method that returns a set of all elements in the molecule
		'''

		return set([a.symbol for a in self.atoms])


	def get_by_element(self, element, blacklist=False):
		'''
		Method that returns a list of atoms corresponding to the atoms in the 
		molecule of a certain element.

		element - string specifying the element

		returns list of indices
		'''

		if blacklist:
			return [a for a in self.atoms if not a.symbol == element]
		return [a for a in self.atoms if a.symbol == element]


	def get_bonded_atoms(self, a, element='any'):
		'''
		Method that returns the indices of atoms bonded to atom a

		a - integer index of atom a

		returns list of integers as indices of atoms bonded to a
		'''

		if element == 'any':
			return a.bonds
		else:
			return [b for b in a.bonds if b.symbol == element]


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
			
			if a.symbol == 'C' and a.hybridisation == 'sp2' and a.HA_valence == 2:
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

		a - atom type
		'''

		self.atoms.remove(a)

		for i in self.atoms:
			try: 
				i.bonds.remove(a)
			except:
				pass

		self.set_bonds()


	def remove_atoms(self, a):
		'''
		Method that removes multiple atoms from the molecule.

		a - atom type
		'''

		for i in a:
			self.remove_atom(i)


	def add_atom(self, element, coords, set_bonds=True):
		'''
		Method that adds an atom to the molecule

		element - string element type of the atom
		coords - np.array with shape (1,3) representing the coordinates of the atom
		'''

		self.atoms.append(Atom(element, coords.flatten()))


	def save(self, file=None, comment='generated via python\n', filetype='xyz'):
		'''
		Method that saves the molecule to a file

		file - string path to new file
		'''

		if file == None:
			file = os.getcwd() + rf'\Molecules\{self.name.lower()}.{filetype}'
		else:
			filetype = file.split('.')[-1]

		if not file.endswith('.' + filetype):
			file += '.' + filetype

		if not comment.endswith('\n'):
			comment += '\n'


		if filetype == 'xyz': #standard xyz format
			with open(file, 'w+') as f:
				f.write(f'{self.natoms}\n')
				f.write(comment)
				for a in self.atoms:
					f.write(f'{a.symbol: <2} \t {a.coords[0]: >8.5f} \t {a.coords[1]: >8.5f} \t {a.coords[2]: >8.5f}\n')

		elif filetype == 'xyzb': #xyz format plus bond information
			with open(file, 'w+') as f:
				f.write(f'{self.natoms}\n')
				f.write(comment)
				for a in self.atoms:
					f.write(f'{a.symbol: <2} \t {a.coords[0]: >8.5f} \t {a.coords[1]: >8.5f} \t {a.coords[2]: >8.5f}\n')

				f.write('\n')
				for a1, a2 in self.get_unique_bonds():
					f.write(f'{self.atoms.index(a1)} \t {self.atoms.index(a2)} \t {a1.bond_orders[a2]}\n')


		utils.message(f'Saved molecule to {file}')

		return file	


	def rotate(self, rotation):
		'''
		Method that rotates atom coordinates by rotation
		'''
		self.rotation += rotation

		r = rotation[0]
		Rx = np.array(([	  1, 	  0,	   0],
					   [	  0, cos(r), -sin(r)],
					   [      0, sin(r),  cos(r)]))

		r = rotation[1]
		Ry = np.array(([ cos(r),  	   0, sin(r)],
					   [ 	  0, 	   1,	   0],
					   [-sin(r), 	   0, cos(r)]))


		for a in self.atoms:
			a.coords = Rx @ Ry @ a.coords

		if hasattr(self, '_elec_stat_pos'):
			self._elec_stat_pos = (Rx @ Ry @ self._elec_stat_pos.T).T

		if hasattr(self, '_dens_pos'):
			for key in self._dens_pos.keys():
				self._dens_pos[key] = (Rx @ Ry @ self._dens_pos[key].T).T