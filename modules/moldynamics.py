import numpy as np
import math
from math import sqrt, exp, pi
import copy
import modules.uff as uff
import scipy.constants as consts
import modules.utils as utils



def boltzmann(x, velocity, mass, temperature):
	return sqrt(mass/(2*pi*consts.k*temperature)) * exp(-mass*velocity**2/(2*consts.k*temperature))


def get_forces(mol, d, ff):
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

	J = np.empty(3*l)

	for i, a in enumerate(atoms):
		for j, c in enumerate(a.coords):
			a.coords[j] += d
			J[3*i+j] = E(mol, morse_potential=True)
			a.coords[j] -= d

	return -(J-e0)/d


def perform_md(mol, ff='uff', run_time=.5e-14, time_step=.5e-15, temperature=273, sample_freq=5):

	if ff == 'uff':
		ff = uff.ForceField()

	atoms = mol.atoms
	l = len(atoms)

	masses = np.expand_dims(np.asarray([a.mass for a in atoms]), 1)

	#initialize velocities
	velocities = np.random.normal(loc=0.5, scale=0.5, size=(l, 3, 12))
	velocities = np.sum(velocities, axis=2) - 6
	# velocities *= np.sqrt(consts.k*temperature/masses)

	print(velocities)

	utils.message(f'Starting molecular dynamics for molecule {mol.name} with {ff.name}.')
	utils.message(f'Run Time: {run_time:.2e} s; Time Step: {time_step:.2e} s; Temperature: {temperature}', 1)

	mol = copy.deepcopy(mol)
	mols = [mol]
	energies = [ff.get_energy(mol)]

	time = 0
	i = 0
	while time < run_time:
		forces = get_forces(mol, 1e-6, ff).reshape((l,3))
		velocities += forces/masses * time_step

		if i%sample_freq == 0:
			mol.center()
			mols.append(copy.deepcopy(mol))
			energies.append(ff.get_energy(mol))

			utils.message((f'Current Time: {time:.2e} with ENERGY = {energies[-1]:.6f} kcal/mol',
						   f'Current Time: {time:.2e} with ENERGY = {energies[-1]:.6f} kcal/mol\nCOORDINATES (angstrom):\n{mol}\n\nVELOCITIES (angstrom/s):\n{velocities}\n'), (1,2))

		for j, a in enumerate(mol.atoms):
			a.coords += velocities[j] * time_step

		time += time_step
		i += 1
	return mols, energies