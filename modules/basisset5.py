import math, json
import numpy as np
import os
from math import pi, exp, sqrt, factorial
from scipy.spatial.distance import cdist


def extended_huckel(molecule, K=1.75):
	'''
	Function that calculates the orbital coefficients via the extended huckel methd
	'''
	atoms = molecule.atoms
	orbitals = []
	for a in atoms:
		for p in a.atomic_orbitals.primitives_valence:
			orbitals.append(p)

	print(orbitals)

	dim = len(orbitals)

	H = np.zeros((dim,dim))

	for i in range(dim):
		for j in range(i, dim):
			if i == j:
				H[i,j] = -orbitals[i].atom.ionisation_energy[sum(orbitals[i].cardinal_powers)] * 0.03676470588235294
			else:
				H[i,j] = H[j,i] = K * overlap_integral(orbitals[i], orbitals[j]) * (H[i,i] + H[j,j])/2

	om = np.zeros((dim,dim))
	for i in range(dim):
		for j in range(i, dim):
			om[i,j] = om[j,i] = overlap_integral(orbitals[i], orbitals[j])


	#construct MO's based on eigenvectors
	mos = []
	energy, weights = np.linalg.eigh(H)
	energy, weights = map(list, zip(*sorted(zip(energy, weights))))

	print(weights[0])
	for e, w in zip(energy, weights):
		mos.append(MolecularOrbital(orbitals, w, e))


	molecule.mos = mos
	return mos


def overlap_matrix(ao):
	om = np.zeros((len(ao), len(ao)))

	for i in range(len(ao)):
		for j in range(i, len(ao)):
			if i == j:
				om[i,j] = 1
			else:
				om[i,j] = om[j,i] = overlap_integral(ao[i], ao[j])


def overlap_integral(ao1, ao2):
	def dfac(n):
		'''
		calculate double factorials n!! = n * (n-2) * ... * (1 or 2)
		'''

		f = 1
		for x in range(n%2, n+1, 2):
			if x != 0:
				f *= x
		return f


	def nCk(n, k):
		#calculates the binomial coefficient
		return factorial(n)/(factorial(k) * factorial(n-k))


	def S(alpha, beta, carda, cardb, A, B):
		s = sqrt(pi/)
		p = alpha + beta
		P = (alpha*A+beta*B)/p
		for i in range(carda):
			for j in range(cardb):
				s += nCk(carda, i) * nCk(cardb, j) * dfac(i+j-1)/((2*p)**((i+j)/2)) * (P-A)**(carda-i) * (P-B)**(cardb-j)


	coorda, coordb = ao1.atom.coords, ao2.atom.coords
	alpha, beta = ao1.a, ao2.a 
	coeffa, coeffb = ao1.c, ao2.c
	carda, cardb = ao1.cardinal_powers, ao2.cardinal_powers

	#loop over all exponents

	s = 0

	for ia, a in enumerate(alpha):
		for ia, b in enumerate(beta):
			p = a+b

			EAB = exp(-a*b/(p) * np.linalg.norm(coorda-coordb)**2)

			Sx = S(a, b, carda[0], cardb[0], coorda[0], coordb[0])
			Sy = S(a, b, carda[1], cardb[1], coorda[1], coordb[1])
			Sz = S(a, b, carda[2], cardb[2], coorda[2], coordb[2])

			s += coeffa[ia] * coeffb[ib] * EAB * Sx * Sy * Sz

	return s



class CGTO:
	'''
	class defining a contracted set of GTO's
	'''

	def __init__(self, centre, a, c, l, m, atom=None):
		self.centre = centre
		self.a = a
		self.l = l
		self.m = m 
		self.c = c
		self.get_cardinal_powers()
		self.norms = [1 for _ in range(len(self.c))]
		self.normalize()
		self.atom = atom
		self.element = atom.symbol
		# print(self.norms)


	def get_cardinal_powers(self):
		if self.l == 0:
			A = (0,0,0)
		if self.l == 1:
			if self.m == -1: A = (1,0,0)
			if self.m ==  0: A = (0,1,0)
			if self.m ==  1: A = (0,0,1)
		if self.l == 2:
			if self.m == -2: A = (2,0,0)
			if self.m == -1: A = (1,1,0)
			if self.m ==  0: A = (1,0,1)
			if self.m ==  1: A = (0,2,0)
			if self.m ==  2: A = (0,1,1)
			if self.m ==  3: A = (0,0,2)
		if self.l == 3:
			if self.m == -3: A = (3,0,0)
			if self.m == -2: A = (2,1,0)
			if self.m == -1: A = (2,0,1)
			if self.m ==  0: A = (1,1,1)
			if self.m ==  1: A = (0,3,0)
			if self.m ==  2: A = (0,2,1)
			if self.m ==  3: A = (0,1,2)
			if self.m ==  4: A = (0,0,3)
		
		self.cardinal_powers = A


	def __repr__(self):
		return f'{self.atom.symbol}({("s", "p", "d", "f")[self.l]}{("x" + str(self.cardinal_powers[0]))*self.cardinal_powers[0]}{("y" + str(self.cardinal_powers[1]))*self.cardinal_powers[1]}{("z" + str(self.cardinal_powers[2]))*self.cardinal_powers[2]})'


	def normalize(self):
		'''
		Method that calculates the norms of the primitives and the contracted GTO set
		'''

		self.norms = []
		for ao in self.primitives_list:
			self.norms.append(sqrt(overlap_integral(ao, ao)))


	def __call__(self, p):
		if len(p.shape) > 1:
			dens = np.zeros(p.shape[0])
		else:
			dens = 0

		l, m, n = self.cardinal_powers
		x, y, z = self.centre
		for a, c, norm in zip(self.a, self.c, self.norms):

			# p = np.expand_dims(np.zeros(p.shape[0]),1)
			px, py, pz = np.hsplit(p, 3)
			hermite = (px-x)**l * (py-y)**m * (pz-z)**n
			r2 = (px-x)**2 + (py-y)**2 + (pz-z)**2
			e = np.maximum(-a * r2, -100)
			dens += (c * hermite * np.exp(e)).flatten()

		return dens


class CGTO:
	def __init__(self, basis_type, atom):
		self.atom = atom
		self.basis_type = basis_type
		self.set_orbitals()


	def set_orbitals(self):
		'''
		Method that returns an orbital function corresponding to n, l, and m
		'''

		self.primitives = dict()
		self.primitives_list = []
		self.primitives_valence = []

		for n in range(len(self.params)):
			self.primitives_valence.append([])
			for l in range(n+1):
				self.primitives[l] = dict()

				for m in ([0], (1,-1,0), (-2,-1,0,1,2,3), (-3,-2,-1,0,1,2,3,4))[l]:
					params = self.params[n]

					exponents = params['exponents']
					try:
						coefficients = params['coefficients'][params['angular_momentum'].index(l)]

						a = [float(a) for a in exponents]
						c = [float(c) for c in coefficients]
						self.primitives[l][m] = CGTO(np.asarray(self.atom.coords), a, c, l, m, self.atom)
						self.primitives_list.append(CGTO(np.asarray(self.atom.coords), a, c, l, m, self.atom))
						self.primitives_valence[-1].append(CGTO(np.asarray(self.atom.coords), a, c, l, m, self.atom))
					
					except:
						raise


		self.primitives_valence = self.primitives_valence[-1]


	@property
	def basis_type(self):
		return self._basis_type


	@basis_type.setter
	def basis_type(self, val):
		'''
		Method that loads basis type whenever the basis type changes.
		'''

		self._basis_type = val
		self.load_basis()


	def load_basis(self):
		'''
		Method that loads the basis set for the atoms in the molecule. If the basis set 
		does not exist in the database, we will download the basis set from https://www.basissetexchange.org/
		using their API
		'''

		bsf_path = os.getcwd()+rf'\Basis_Sets\{self.basis_type}.bsf'
		if not os.path.exists(bsf_path):
			print(f'Error: Basis set {self.basis_type} not found, downloading ...')
			import requests
			response = requests.get("http://basissetexchange.org" + f'/api/basis/{self.basis_type}/format/json')
			if response:
				print('Succesfully obtained basis set file')
				with open(bsf_path, 'w+') as f:
					f.write(response.text)
				self.load_basis()
			else:
				print('Failed to obtain basis set file')
		else:
			print(f'Succesfully loaded {self.basis_type}.bsf')
			with open(bsf_path, 'r') as f:
				self.params = json.load(f)['elements'][str(self.atom.atomic_number)]['electron_shells']


class MolecularOrbital:
	def __init__(self, atomic_orbitals, weights, energy=None):
		self.atomic_orbitals = atomic_orbitals
		self.weights = weights
		self.energy = energy


	def __call__(self, p):
		dens = np.zeros(p.shape[0])
		for ao, w in zip(self.atomic_orbitals, self.weights):
			dens += ao(p) * w

		return dens



class BasisLoader:
	def __init__(self, molecule, basis_type):
		self.molecule = molecule
		self.basis_type = basis_type

	@property
	def basis_type(self):
		return self._basis_type


	@basis_type.setter
	def basis_type(self, val):
		'''
		Method that loads basis type whenever the basis type changes.
		'''

		self._basis_type = val
		self.load_basis()


	def load_basis(self):
		'''
		Method that loads the basis set for the atoms in the molecule. If the basis set 
		does not exist in the database, we will download the basis set from https://www.basissetexchange.org/
		using their API
		'''

		bsf_path = os.getcwd()+rf'\Basis_Sets\{self.basis_type}.bsf'
		if not os.path.exists(bsf_path):
			print(f'Error: Basis set {self.basis_type} not found, downloading ...')
			import requests
			response = requests.get("http://basissetexchange.org" + f'/api/basis/{self.basis_type}/format/json')
			if response:
				print('Succesfully obtained basis set file')
				with open(bsf_path, 'w+') as f:
					f.write(response.text)
				self.load_basis()
			else:
				print('Failed to obtain basis set file')
		else:
			print(f'Succesfully loaded {self.basis_type}')
			with open(bsf_path, 'r') as f:
				# self.params = json.load(f)['elements'][str(self.atom.atomic_number)]['electron_shells']
				self.params = json.load(f)['elements']


	def set_orbitals(self):
		for atom in self.molecule.atoms:

			atom.primitives = {}
			atom.primitives_list = []
			atom.primitives_valence = []

			atom_params = self.params[str(atom.atomic_number)]['electron_shells']
			for n in range(len(atom_params)):
				atom.primitives_valence.append([])
				for l in range(n+1):
					atom.primitives[l] = dict()

					for m in ([0], (-1,0,1), (-2,-1,0,1,2,3), (-3,-2,-1,0,1,2,3,4))[l]:
						params = atom_params[n]
						exponents = params['exponents']

						try:
							coefficients = params['coefficients'][params['angular_momentum'].index(l)]

							a = [float(a) for a in exponents]
							c = [float(c) for c in coefficients]
							atom.primitives[l][m] = CGTO(np.asarray(atom.coords), a, c, l, m, atom)
							atom.primitives_list.append(CGTO(np.asarray(atom.coords), a, c, l, m, atom))
							atom.primitives_valence[-1].append(CGTO(np.asarray(atom.coords), a, c, l, m, atom))
						
						except:
							raise

		self.primitives_valence = self.primitives_valence[-1]