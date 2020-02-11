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


def minimize(mol, ff=uff.ForceField(), max_steps=1500, converge_thresh=5e-2, step_factor=4e-4):
	utils.message(f'Starting geometry optimisation for molecule {mol.name} using {ff.name}.')
	utils.message(f'Max. Steps: {max_steps}; Step-Factor: {step_factor:.2e}; Converge thresh.: {converge_thresh:.2e}', 1)

	mol = copy.deepcopy(mol)
	mols = [mol]
	
	for i in range(max_steps):
		forces = get_forces(mol, 1e-6, ff)
		if np.all(np.absolute(forces) < converge_thresh):
			break

		for j, a in enumerate(mol.atoms):
			a.coords += forces[j] * step_factor

		if i%10 == 0:
			mol.center()
			mols.append(copy.deepcopy(mol))

			utils.message(f'Current Step: {i} with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n', 2)

	if i < max_steps-1:
		utils.message(f'Molecule optimization succesful after {i+1} steps.', 1, colour='green')
		utils.message(f'ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n', 2, colour='green')
	else:
		utils.message(f'Molecule optimization failed after {i+1} steps.', 1, colour='red')
		utils.message(f'ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n', 2, colour='red')
	
	return mols