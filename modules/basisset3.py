import numpy as np
import os
import gbasis.parsers as gb_pars



class BasisSet:
	def __init__(self, atoms, basis_type='STO-6G'):
		self.atoms = atoms
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

		gbs_path = os.getcwd()+rf'\Basis_Sets\{self.basis_type}.gbs'
		if not os.path.exists(gbs_path):
			print(f'Error: Basis set {self.basis_type} not found, downloading ...')
			import requests
			print(requests.get("http://basissetexchange.org/api/formats").text)
			response = requests.get("http://basissetexchange.org" + f'/api/basis/{self.basis_type}/format/Gaussian94')
			if response:
				print('Succesfully obtained basis set file')
				with open(gbs_path, 'w+') as f:
					f.write(response.text)
				self.load_basis()
			else:
				print('Failed to obtain basis set file')
		else:
			print(f'Succesfully loaded {self.basis_type}.gbs')
			all_basis_dict = gb_pars.parse_gbs(gbs_path)
			self.basis = gb_pars.make_contractions(all_basis_dict, 
						[a.symbol for a in self.atoms], 
						np.asarray([a.coords for a in self.atoms]))


