import numpy as np
import math

class Minimizer:
	def __init__(self, forcefield, molecule):
		self.ff = forcefield
		self.mol = molecule


	#maths
	def get_jacobian(self):
		'''
		Method that returns the jacobian vector of the molecule.
		The jacobian vector is a vector containing the derivatives
		of the energy function with respect to some coordinates.
		For simplicity we shall take the positions of the atoms
		as coordinates. Bond lengths and bond angles would be 
		better, but also more complex.
		'''

		E = self.ff.get_energy
		m = self.mol

		d = 10**-8
		J1 = np.empty(len(m.atoms)*3)
		J2 = np.empty(len(m.atoms)*3)
		for i, a in enumerate(m.atoms):
			for j, c in enumerate(a.coords):
				e0 = E(m)
				a.coords[j] += d
				e1 = E(m)
				a.coords[j] -= 2*d
				e2 = E(m)
				a.coords[j] += d
				
				J1[i+j] = (e1-e0)/d
				J2[i+j] = (e2-e0)/d

		return J1

	def minimize(self):
		for i in range(10):
			J = self.get_jacobian()
			print(J)
			for i, a in enumerate(self.mol.atoms):
				a.coords -= J[i*3:i*3+3] * 0.01

			print(self.mol.atoms[0])
