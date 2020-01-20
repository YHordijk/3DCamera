import numpy as np
import os
import openfermion as of
import openfermion.hamiltonians as ofh





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

		geo = list(zip([a.symbol for a in self.atoms], [a.coords for a in self.atoms]))
		print(geo)
		self.basis = ofh.MolecularData(geometry=geo, basis=self.basis_type, multiplicity=1)


