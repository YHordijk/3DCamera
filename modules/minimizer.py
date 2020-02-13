import numpy as np
import math
import copy
import modules.uff as uff
import modules.utils as utils


def get_forces(mol, d, ff, use_torsions):
	'''
	Method that calculates the forces acting
	on the atoms in mol. Uses atomic 
	cartesian coordinates to calculate the 
	forces.

	Returns a nx3 matrix containing the forces.
	'''

	E = ff.get_energy
	e0 = E(mol)

	atoms = mol.atoms
	l = len(atoms)

	if use_torsions:
		torsions = list(mol.get_rotatable_bonds())

		# J = np.empty(len(torsions))

		J = np.empty(3*l + len(torsions))

		for i, a in enumerate(atoms):
			for j, c in enumerate(a.coords):
				a.coords[j] += d
				J[3*i+j] = E(mol, morse_potential=True)
				a.coords[j] -= d

		for k, t in enumerate(torsions):
			mol.rotate_bond(t[0],t[1], d)
			J[3*l + k] = E(mol, morse_potential=True)

	else:
		J = np.empty(3*l)

		for i, a in enumerate(atoms):
			for j, c in enumerate(a.coords):
				a.coords[j] += d
				J[3*i+j] = E(mol, morse_potential=True)
				a.coords[j] -= d

	return -(J-e0)/d


def minimize(mol, ff='uff', max_steps=1500, converge_thresh=8e-2, step_factor=4e-4, sample_freq=10, use_torsions=True, max_step_size=0.3, method='sd'):
	'''
	Energy minimization method that attempts to optimize the structureof a molecule using a given force field
	(must have a get_energy() method). The method used is the steepest descent method, which follows the 
	negative energy gradient with respect to the coordinates of the molecule multiplied by step_factor.
	Every sample_freq iterations a snapshot is taken.
	'''

	assert(max_steps > 0)
	assert(converge_thresh > 0)
	assert(step_factor > 0)
	assert(sample_freq > 0)

	if ff == 'uff':
		ff = uff.ForceField()

	utils.message(f'Starting geometry optimisation for molecule {mol.name} with {ff.name} using steepest descent.')
	utils.message(f'Max. Steps: {max_steps}; Step-Factor: {step_factor:.2e}; Max. Step Size: {max_step_size}; Converge Thresh.: {converge_thresh:.2e}; Apply to Torsion: {use_torsions}', 1)

	mol = copy.deepcopy(mol)
	mols = [mol]
	print(mols[0])
	energies = [ff.get_energy(mol)]

	if method == 'sd':
		for i in range(max_steps):
			forces = get_forces(mol, 1e-6, ff, use_torsions=use_torsions)

			if i%sample_freq == 0 or np.all(np.absolute(forces) < converge_thresh):
				mol.center()
				mols.append(copy.deepcopy(mol))
				energies.append(ff.get_energy(mol))
				if np.all(np.absolute(forces) < converge_thresh):
					break

				utils.message((f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol',
							   f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (1,2))

			if use_torsions:
				for j, a in enumerate(mol.atoms):
					a.coords += (forces[3*j:3*j+3] * step_factor)
				for k, t in enumerate(list(mol.get_rotatable_bonds())):
					mol.rotate_bond(t[0],t[1], step_factor * forces[3*len(mol.atoms) + k])
			else:
				for j, a in enumerate(mol.atoms):
					a.coords += step_sizes[3*j:3*j+3]

			
		if i < max_steps-1:
			utils.message((f'Molecule optimization succesful after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol',
						   f'Molecule optimization succesful after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (0,2), colour='green')
		else:
			utils.message((f'Molecule optimization failed after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol',
						   f'Molecule optimization failed after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (0,2), colour='red')
			
		return mols, energies

	if method == 'cg':
		#first step is the same as steepest descent:
		forces = get_forces(mol, 1e-6, ff, use_torsions=use_torsions)
		step = forces * step_factor
		prev_step = step

		if use_torsions:
			for j, a in enumerate(mol.atoms):
				a.coords += (step[3*j:3*j+3])
			for k, t in enumerate(list(mol.get_rotatable_bonds())):
				mol.rotate_bond(t[0],t[1], step[3*len(mol.atoms) + k])
		else:
			for j, a in enumerate(mol.atoms):
				a.coords += step[3*j:3*j+3]

		#next steps:
		for i in range(max_steps):
			forces = get_forces(mol, 1e-6, ff, use_torsions=use_torsions)
			gamma = np.dot(-forces, -forces)/np.dot(prev_step, prev_step)
			step = forces * gamma*prev_step
			prev_step = step

			if i%sample_freq == 0 or np.all(np.absolute(forces) < converge_thresh):
				mol.center()
				mols.append(copy.deepcopy(mol))
				energies.append(ff.get_energy(mol))
				if np.all(np.absolute(forces) < converge_thresh):
					break

				utils.message((f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol',
							   f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (1,2))


			if use_torsions:
				for j, a in enumerate(mol.atoms):
					a.coords += (step[3*j:3*j+3])
				for k, t in enumerate(list(mol.get_rotatable_bonds())):
					mol.rotate_bond(t[0],t[1], step[3*len(mol.atoms) + k])
			else:
				for j, a in enumerate(mol.atoms):
					a.coords += step[3*j:3*j+3]

			return mols, energies
