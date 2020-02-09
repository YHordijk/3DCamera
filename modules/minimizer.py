import numpy as np
import math
import copy
import modules.uff as uff

#maths
def get_forces(mol, d=1e-7, ff=uff.ForceField()):
	'''
	Method that calculates the forces acting
	on the atoms in mol. Uses atomic 
	cartesian coordinates to calculate the 
	forces.

	Returns a nx3 matrix containing the forces.
	'''

	E = ff.get_energy
	e0 = E(mol)

	J = np.empty((len(mol.atoms), 3))

	for i, a in enumerate(mol.atoms):
		for j, c in enumerate(a.coords):
			a.coords[j] += d
			J[i,j] = E(mol)
			a.coords[j] -= d

	return -(J-e0)/d


def minimize(mol, ff=uff.ForceField(), steps=100, converge=1e-3, step_factor=1e-3):
	for i in range(steps):
		forces = get_forces(mol, 1e-7, ff)
		for j, a in enumerate(mol.atoms):
			a.coords += forces[j] * step_factor
		print(i/steps*100)
		print(a.coords)

	return mol