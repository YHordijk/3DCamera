import numpy as np
import math
import copy
import modules.uff as uff
import modules.utils as utils

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
			J[i,j] = E(mol, morse_potential=True)
			a.coords[j] -= d

	return -(J-e0)/d


def minimize(mol, ff=uff.ForceField(), steps=1500, converge=1e-3, step_factor=1e-4, copy_mol=True):
	utils.message(f'Starting geometry optimisation for molecule {mol.name} using {ff.name}.')
	utils.message(f'Steps: {steps}, Step-Factor: {step_factor}', 1)
	if copy_mol: mol = copy.deepcopy(mol)
	for i in range(steps):
		if i%100 == 0:
			utils.message(f'Progress: {int(i/steps*100)}%', 2)
		forces = get_forces(mol, 1e-6, ff)
		for j, a in enumerate(mol.atoms):
			a.coords += forces[j] * step_factor
	utils.message('Progress: 100%', 2)
	return mol