import math, json
import numpy as np
import os
from math import pi



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


	def __getitem__(self, val):
		'''
		Support for calling the object, returns the parameters.
		'''

		return self.params[val]


	def orbital_func(self, l, ml, a, c, x, y, z, r2):
		e = np.exp(-a*r2)
		e[np.abs(e) < 0.0001] = 0
		d = c * (2*a/pi)**(3/4) * e

		if l == 0:
			return d
		if l == 1:
			if ml == -1: return lambda x, y, z: d * x
			if ml ==  0: return lambda x, y, z: d * y
			if ml ==  1: return lambda x, y, z: d * z
		if l == 2:
			if ml == -2: return lambda x, y, z: d * (x**2-y**2)
			if ml == -1: return lambda x, y, z: d * x*y
			if ml ==  0: return lambda x, y, z: d * x*z
			if ml ==  1: return lambda x, y, z: d * y*z
			if ml ==  2: return lambda x, y, z: d * z*z
		if l == 3:
			if ml == -3: return lambda x, y, z: d * y*(y*y-3*x*x)
			if ml == -2: return lambda x, y, z: d * x*(x*x-3*y*y)
			if ml == -1: return lambda x, y, z: d * z*(x*x-y*y)
			if ml ==  0: return lambda x, y, z: d * x*y*z
			if ml ==  1: return lambda x, y, z: d * y*z*z
			if ml ==  2: return lambda x, y, z: d * x*z*z
			if ml ==  3: return lambda x, y, z: d * z*z*z
		return 0


	def orbital(self, l, ml, a, c, x, y, z, r2):
		e = np.exp(-a*r2)
		e[np.abs(e) < 0.0001] = 0
		d = c * (2*a/pi)**(3/4) * e

		if l == 0:
			return d
		if l == 1:
			if ml == -1: return d * x
			if ml ==  0: return d * y
			if ml ==  1: return d * z
		if l == 2:
			if ml == -2: return d * (x**2-y**2)
			if ml == -1: return d * x*y
			if ml ==  0: return d * x*z
			if ml ==  1: return d * y*z
			if ml ==  2: return d * z*z
		if l == 3:
			if ml == -3: return d * y*(y*y-3*x*x)
			if ml == -2: return d * x*(x*x-3*y*y)
			if ml == -1: return d * z*(x*x-y*y)
			if ml ==  0: return d * x*y*z
			if ml ==  1: return d * y*z*z
			if ml ==  2: return d * x*z*z
			if ml ==  3: return d * z*z*z
		return 0



	def __call__(self, p, mo=0):
		'''
		Support for calling the object like a function, takes position as argument
		Returns density of basis set at the location
		'''

		d = 0
		for i, atom in enumerate(self.atoms):
			x, y, z = np.hsplit(p - atom.coords, 3)
			r2 = x**2 + y**2 + z**2
			#loop over all atoms and get the correct basis set parameters depending on the element
			params = self.params[str(atom.atomic_number)]['electron_shells']
			# for param in params:
			# 	#find the last (valence) parameters of the element
			# 	params = param
			params = params[-1]
			#get the exponential and coefficient parameters
			expon = params['exponents']
			coeff = params['coefficients']
			angmom = params['angular_momentum']

			#implement formula d = Ne^(-ar^2), N = (2a/pi)^(3/4)

			for l, coef in zip(angmom, coeff):
				for a, c in zip(expon, coef):
					# return self.orbital(3,-1, a, c, x, y, z, r2)

					if mo == 1:
						d += self.orbital(0,0, float(a), float(c), x, y, z, r2)
					if mo == 8:
						if i == 0:
							d += self.orbital(0, 0, float(a), float(c), x, y, z, r2)
						else:
							d -= self.orbital(0, 0, float(a), float(c), x, y, z, r2) 
					if mo == 6:
						if i == 0:
							d -= self.orbital(1, 1, float(a), float(c), x, y, z, r2)
						else:
							d += self.orbital(0, 0, float(a), float(c), x, y, z, r2) * (-1,1)[i in (1,3)]
					else:
						d += sum([self.orbital(l, ml, float(a), float(c), x, y, z, r2) for ml in range(-l, l+1)])
					
		return d**2


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
				self.params = json.load(f)['elements']



	
