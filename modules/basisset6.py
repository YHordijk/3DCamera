import math, json
import numpy as np
import os
from math import pi, exp, sqrt, factorial
from scipy.spatial.distance import euclidean



def extended_huckel(molecule, K=1.75):
	'''
	Function that performs the extended huckel method

	Returns energies and mo's
	'''

	aos = molecule.basis.atomic_orbitals

	dim = len(aos)
	H = np.zeros((dim, dim))

	for i in range(dim):
		H[i,i] = -aos[i].atom.ionisation_energy[sum(aos[i].cardinality)] * 0.036749405469679

	for i in range(dim):
		for j in range(i+1, dim):
			H[i,j] = H[j,i] = K * overlap_integral(aos[i], aos[j]) * (H[i,i] + H[j,j])/2


	energies, weights = np.linalg.eigh(H)

	molecule.molecular_orbitals = []
	for energy, weight in zip(energies, weights.T):
		molecule.molecular_orbitals.append(MolecularOrbital(aos, weight, energy))


def overlap_matrix(molecule):
	basis = molecule.basis
	ao = basis.atomic_orbitals
	om = np.zeros((len(ao), len(ao)))

	for i in range(len(ao)):
		for j in range(i, len(ao)):
			om[i,j] = om[j,i] = overlap_integral(ao[i], ao[j])

	return om


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
		
		p = alpha + beta
		P = (alpha*A+beta*B)/p

		s_pre = sqrt(pi/p)
		s = 0
		for i in range(carda+1):
			for j in range(cardb+1):
				if (i+j)%2 == 0:
					s += nCk(carda, i) * nCk(cardb, j) * dfac(i+j-1)/((2*p)**((i+j)/2)) * (P-A)**(carda-i) * (P-B)**(cardb-j)

		return s_pre * s

	coorda, coordb = ao1.centre, ao2.centre
	alpha, beta = ao1.exponents, ao2.exponents
	coeffa, coeffb = ao1.coefficients, ao2.coefficients
	carda, cardb = ao1.cardinality, ao2.cardinality

	#loop over all exponents

	s = 0

	for ia, a in enumerate(alpha):
		for ib, b in enumerate(beta):
			p = a+b

			EAB = exp(-a*b/(p) * np.linalg.norm(coorda-coordb)**2)

			Sx = S(a, b, carda[0], cardb[0], coorda[0], coordb[0])
			Sy = S(a, b, carda[1], cardb[1], coorda[1], coordb[1])
			Sz = S(a, b, carda[2], cardb[2], coorda[2], coordb[2])

			s += coeffa[ia] * coeffb[ib] * EAB * Sx * Sy * Sz * ao1.norm * ao2.norm

	return s



class Basis:
	'''
	Master class containing all information for basis set of a molecule
	'''

	def __init__(self, molecule, basis_type):
		self.molecule = molecule
		self.basis_type = basis_type
		self.set_orbitals()

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


	def get_cardinal_powers(self, l, m):
		if l == 0:
			return (0,0,0)
		if l == 1:
			if m == -1: return (1,0,0)
			if m ==  0: return (0,1,0)
			if m ==  1: return (0,0,1)
		if l == 2:
			if m == -2: return (2,0,0)
			if m == -1: return (1,1,0)
			if m ==  0: return (1,0,1)
			if m ==  1: return (0,2,0)
			if m ==  2: return (0,1,1)
			if m ==  3: return (0,0,2)
		if l == 3:
			if m == -3: return (3,0,0)
			if m == -2: return (2,1,0)
			if m == -1: return (2,0,1)
			if m ==  0: return (1,1,1)
			if m ==  1: return (0,3,0)
			if m ==  2: return (0,2,1)
			if m ==  3: return (0,1,2)
			if m ==  4: return (0,0,3)
		

	def set_orbitals(self):
		self.atomic_orbitals = []
		for atom in self.molecule.atoms:
			atom_params = self.params[str(atom.atomic_number)]['electron_shells']

			for n in range(len(atom_params)):
				for l in range(n+1):
					for m in ([0], (1, -1,0), (-2,-1,0,1,2,3), (-3,-2,-1,0,1,2,3,4))[l]:
						params = atom_params[n]
						exponents = params['exponents']

						try:
							coefficients = params['coefficients'][params['angular_momentum'].index(l)]

							a = [float(a) for a in exponents]
							c = [float(c) for c in coefficients]

							if n == len(atom_params) -1:
								self.atomic_orbitals.append(AtomicOrbital(atom.coords, n+1, a, c, self.get_cardinal_powers(l, m), atom))

						except:
							pass


class AtomicOrbital:
	'''
	Class containing the info for atomic orbitals for an atom
	'''

	def __init__(self, centre, principal_qn, exponents, coefficients, cardinality, atom=None):
		self.centre = np.asarray(centre)
		self.principal = principal_qn
		self.exponents = exponents
		self.coefficients = coefficients
		self.cardinality = cardinality
		self.atom = atom

		self.norm = 1
		self.normalize()


	def __repr__(self):
		l = sum(self.cardinality)
		c = self.cardinality
		return f'{self.atom.symbol}({self.principal}{("s", "p", "d", "f")[l]}{("x" + str(c[0])*(c[0]>1))*c[0]}{("y" + str(c[1])*(c[1]>1))*c[1]}{("z" + str(c[2])*(c[2]>1))*c[2]})'


	def normalize(self):
		self.norm = 1/sqrt(overlap_integral(self, self))


	def evaluate(self, p):
		Ax, Ay, Az = self.centre
		x, y, z = np.hsplit(p.T,3)
		x, y, z = x.flatten(), y.flatten(), z.flatten()
		ax, ay, az = self.cardinality

		dens = np.zeros(p.size//3)

		for a, c in zip(self.exponents, self.coefficients):
			dens += c * (x-Ax)**ax * (y-Ay)**ay * (z-Az)**az * np.exp(-a*np.linalg.norm(p.T - self.centre, axis=1)**2) 
		return dens* self.norm



class MolecularOrbital:
	def __init__(self, aos, weights, energy):
		self.aos = aos 
		self.weights = weights
		self.energy = energy

	def evaluate(self, p):
		dens = np.zeros(p.size//3)
		for ao, w in zip(self.aos, self.weights):
			dens += w * ao.evaluate(p)

		return dens

