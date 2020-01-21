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

	print(om)

	# #calculate huckelmatrix
	# dim = om.shape[0]
	# hm = np.zeros((dim,dim))

	# #first set diagonals elements
	# for i in range(dim):
	# 	hm[i,i] = -orbitals[i].atom.ionisation_energy[sum(orbitals[i].cardinal_powers)]

	# print(om)
	# print()
	# print(hm)

	# #set the off-diagonals
	# for i in range(1,dim):
	# 	for j in range(i+1, dim):
	# 		#H[i,j] = K * S[i,j] * (H[i,i] + H[j,j])/2
	# 		#H[i,i] is atomisation energy corresponding to either s or p orbitals

	# 		hm[i,j] = hm[j,i] = -K * om[i,j] * (hm[i,i] + hm[j,j])/2
	# 		if i==1 and j==4:
	# 			print()
	# 			print(-K * om[i,j] * (hm[i,i] + hm[j,j])/2)





	print(H)
	#construct MO's based on eigenvectors
	mos = []
	energy, weights = np.linalg.eigh(H)
	energy, weights = map(list, zip(*sorted(zip(energy, weights))))

	print(weights[0])
	for e, w in zip(energy, weights):
		mos.append(MolecularOrbital(orbitals, w, e))


	molecule.mos = mos
	return mos





def overlap_integral(ao1, ao2):
	'''
	Function that calculates the overlap between two gaussians
	based on https://joshuagoings.com/2017/04/28/integrals/
	'''

	def dfac(n):
		'''
		calculate double factorials n!! = n * (n-2) * ... * (1 or 2)
		'''

		f = 1
		for x in range(n%2, n+1, 2):
			if x != 0:
				f *= x
		return f


	def E(i, j, t, Qx, a, b):
		p = a + b
		q = a*b/p

		if t < 0 or t > i+j:
			return 0
		elif i == j == t == 0:
			return exp(-q*Qx**2)
		elif i == 0:
			return 1/(2*p) * E(i, j-1, t-1, Qx, a, b) + q*Qx/b * E(i, j-1, t, Qx, a, b) + (t+1) * E(i, j-1, t+1, Qx, a, b)
		else:
			return 1/(2*p) * E(i-1, j, t-1, Qx, a, b) + q*Qx/b * E(i-1, j, t, Qx, a, b) + (t+1) * E(i-1, j, t+1, Qx, a, b)

	def overlap(a, A, lmn1, b, B, lmn2):
		l1, m1, n1 = lmn1
		l2, m2, n2 = lmn2

		S1 = E(l1, l2, 0, A[0]-B[0], a, b)
		S2 = E(m1, m2, 0, A[1]-B[1], a, b)
		S3 = E(n1, n2, 0, A[2]-B[2], a, b)

		return S1 * S2 * S3 * (pi/(a+b))**1.5


	#calculate the integral:
	s = 0
	for ia, ca in enumerate(ao1.c):
		for ib, cb in enumerate(ao2.c):
			s += ao1.norms[ia] * ao2.norms[ib] * ca * cb * \
				overlap(ao1.a[ia], ao1.centre, ao1.cardinal_powers,
						ao2.a[ib], ao2.centre, ao2.cardinal_powers)


	return s


class BasisSet:
	def __init__(self, basis_type, atoms=[]):
		self.basis_type = basis_type
		self.atoms = atoms


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

		def dfac(n):
			'''
			calculate double factorials n!! = n * (n-2) * ... * (1 or 2)
			'''
			f = 1
			for x in range(n%2, n+1, 2):
				if x != 0:
					f *= x
			return f


		l, m, n = self.cardinal_powers
		L = l+m+n

		# self.norms = np.sqrt(np.power(2,2*(L)+1.5)*
		# 	np.power(self.a,L+1.5)/
		# 	dfac(2*l-1)/dfac(2*m-1)/
		# 	dfac(2*n-1)/np.power(np.pi,1.5))

		self.norms = []

		for a in self.a:
			self.norms.append(np.sqrt(np.sqrt(2*a/pi) * 2**(2*L) * a**L / \
				(dfac(2*l-1) * dfac(2*m-1) * dfac(2*n-1))))

		prefactor = np.power(np.pi,1.5)* \
			dfac(2*l - 1)*dfac(2*m - 1)*dfac(2*n - 1)/np.power(2.0,L)


		N = 0.0
		num_exps = len(self.a)
		for ia in range(num_exps):
			for ib in range(num_exps):
				N += self.norms[ia]*self.norms[ib]*self.c[ia]*self.c[ib]/\
						 np.power(self.a[ia] + self.a[ib],L+1.5)

		N *= prefactor
		N = np.power(N,-0.5)
		for ia in range(num_exps):
			self.c[ia] *= N


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


class AtomicOrbital(BasisSet):
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