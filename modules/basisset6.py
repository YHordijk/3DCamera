import math, json
import numpy as np
import os
from math import pi, exp, sqrt, factorial
from scipy.spatial.distance import cdist



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