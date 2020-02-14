import numpy as np
import math
import copy
import modules.uff as uff
import modules.utils as utils
import modules.molecule6 as mol6
import numdifftools as nd


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

			if (i+1)%sample_freq == 0 or np.all(np.absolute(forces) < converge_thresh):
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
						   f'Molecule optimization succesful after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}'), (0,2), colour='green')
		else:
			utils.message((f'Molecule optimization failed after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol',
						   f'Molecule optimization failed after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}'), (0,2), colour='red')
			
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

			if (i+1)%sample_freq == 0 or np.all(np.absolute(forces) < converge_thresh):
				mol.center()
				mols.append(copy.deepcopy(mol))
				energies.append(ff.get_energy(mol))
				if np.all(np.absolute(forces) < converge_thresh):
					break

				utils.message((f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol',
							   f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}'), (1,2))


			if use_torsions:
				for j, a in enumerate(mol.atoms):
					a.coords += (step[3*j:3*j+3])
				for k, t in enumerate(list(mol.get_rotatable_bonds())):
					mol.rotate_bond(t[0],t[1], step[3*len(mol.atoms) + k])
			else:
				for j, a in enumerate(mol.atoms):
					a.coords += step[3*j:3*j+3]

			return mols, energies


class Minimizer:

	def __init__(self, mol, ff='uff'):
		self._mol = mol
		self._n_torsions = len(list(mol.get_rotatable_bonds()))
		self._n_atoms = len(mol.atoms)

		if ff == 'uff':
			self._ff = uff.ForceField()

	def _apply_coords(self, coords, use_torsions=False):
		'''
		Method that handles setting the coords of the molecule.
		'''

		n = self._n_atoms

		#update cartesian coords
		for c, a in zip(coords[:n*3].reshape((n,3)), self._mol.atoms):
			a.coords = c

		#update torsions
		if use_torsions: 
			for t, a in zip(coords[n*3:], self._mol.get_rotatable_bonds()):
				self.mol.rotate_bond(*a, t)


	def _get_energy_from_coords(self, coords):
		self._apply_coords(coords)
		return self._ff.get_energy(self._mol)



	def minimize(self, method='sd', max_steps=1500, step_factor=4e-4, sample_freq=10, use_torsions=True, converge_thresh=8e-2):
		'''
		Method that performs the energy minimization of self.mol

		method: str - specify method (steepest descent 'sd', conjugate gradient 'cg')
		max_steps: int - specify the maximum number of iterations allowed for minimization
		sample_freq: int - specify the frequency with which a snapshot is taken (in iterations)
		use_torsions: bool - specify whether to use rotation around bonds as coordinates in minimization

		returns tuple of lists containing molecule objects and energies
		'''

		assert(max_steps > 0)
		assert(step_factor > 0)
		assert(sample_freq > 0)

		mol = self._mol

		utils.message((f'Starting geometry optimisation for molecule {mol.name} with {self._ff.name} using steepest descent.',
					   f'Starting geometry optimisation for molecule {mol.name} with {self._ff.name} using steepest descent.\nINITIAL COORDINATES (angstrom):\n{mol}'), (0,1))
		utils.message(f'Max. Steps: {max_steps}; Step-Factor: {step_factor:.2e}; Convergence Thresh.: {converge_thresh:.2e}; Apply to Bond Rotation: {use_torsions}', 1)


		energy = self._ff.get_energy
		mols = [copy.deepcopy(mol)]
		energies = [energy(mol)]


		if method == 'sd':
			gradient = nd.Gradient(self._get_energy_from_coords)


			#get the coordinates:
			coords = [a.coords for a in mol.atoms]  #cartesian coordinates of the atoms
			if use_torsions:
				coords += [0 for _ in mol.get_rotatable_bonds()] #torsions
			coords = np.asarray(coords).flatten().astype(float) #flatten vector

			

			#start minimization
			for i in range(max_steps):
				#calculate gradients and change coords
				forces = -gradient(coords) * step_factor
				coords += forces

				converged = np.all(np.absolute(forces) < converge_thresh)

				#sample mol
				if (i+1)%sample_freq == 0 or converged:

					mols.append(copy.deepcopy(mol))
					energies.append(energy(mol))

					if converged: break

					utils.message((f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol',
								   f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol\nCOORDINATES (angstrom):\n{mol}\n\nFORCES (kcal/mol/angstrom):\n{forces}\n'), (1,2))

			#done
			if i < max_steps-1:
				utils.message((f'Molecule optimization succesful after {i+1} steps with ENERGY = {energy(mol):.6f} kcal/mol',
							   f'Molecule optimization succesful after {i+1} steps with ENERGY = {energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (0,2), colour='green')
			else:
				utils.message((f'Molecule optimization failed after {i+1} steps with ENERGY = {energy(mol):.6f} kcal/mol',
							   f'Molecule optimization failed after {i+1} steps with ENERGY = {energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (0,2), colour='red')
				
			return mols, energies



# def minimize2(mol, ff='uff', max_steps=1500, converge_thresh=8e-2, step_factor=4e-4, sample_freq=10, use_torsions=True, method='sd'):
# 	'''
# 	doc-string
# 	'''

# 	assert(max_steps > 0)
# 	assert(converge_thresh > 0)
# 	assert(step_factor > 0)
# 	assert(sample_freq > 0)

# 	if ff == 'uff':
# 		ff = uff.ForceField()

# 	utils.message((f'Starting geometry optimisation for molecule {mol.name} with {ff.name} using steepest descent.',
# 				   f'Starting geometry optimisation for molecule {mol.name} with {ff.name} using steepest descent.\nINITIAL COORDINATES (angstrom):\n{mol}'), (0,1))
# 	utils.message(f'Max. Steps: {max_steps}; Step-Factor: {step_factor:.2e}; Max. Step Size: {max_step_size}; Converge Thresh.: {converge_thresh:.2e}; Apply to Torsion: {use_torsions}', 1)

# 	elements = [a.symbol for a in mol.atoms]
# 	mols = [copy.deepcopy(mol)]

# 	if method == 'sd':
# 		gradient = nd.Gradient(ff.get_energy)

# 		#get the coordinates:
# 		coords = [a.coords for a in mol.atoms]  #cartesian coordinates of the atoms
# 		if use_torsions:
# 			coords += [mol.get_rotatable_bonds] #torsions

# 		coords = np.asarray(coords).flatten() #make vector

# 		for i in range(max_steps):
# 			#calculate gradients and change coords
# 			forces = -gradient(coords) * step_factor
# 			coords += forces

# 			if (i+1)%sample_freq == 0 or converged:
# 				for a in mol.atoms:
# 					a.coords = coords[3*j:3*j+3]
# 				mols.append(copy.deepcopy(mol))
# 				energies.append(ff.get_energy(mol))

# 				if converged: break

# 				utils.message((f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol',
# 							   f'Current Step: {i+1} with ENERGY = {energies[-1]:.6f} kcal/mol\nCOORDINATES (angstrom):\n{mol}\n\nFORCES (kcal/mol/angstrom):\n{forces}\n'), (1,2))

# 		if i < max_steps-1:
# 			utils.message((f'Molecule optimization succesful after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol',
# 						   f'Molecule optimization succesful after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (0,2), colour='green')
# 		else:
# 			utils.message((f'Molecule optimization failed after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol',
# 						   f'Molecule optimization failed after {i+1} steps with ENERGY = {ff.get_energy(mol):.6f} kcal/mol\nCOORDINATES (angstrom):\n\n{mol}\n\nFORCES (kcal/mol/angstrom):\n\n{forces}\n'), (0,2), colour='red')
			
# 		return mols, energies
