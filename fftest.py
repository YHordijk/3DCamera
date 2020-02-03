import modules.molecule6 as mol 
import modules.reaxff as reaxff
import modules.uff as uff
import modules.utils as utils
import os
import numpy as np

utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(True)


ff = uff.ForceField()
ethane = mol.Molecule('ethane')
ethane_coords = np.asarray([a.coords for a in ethane.atoms])

for r in np.arange(1.3, 1.6, 0.01):
	coords = ethane_coords.copy()
	coords -= np.array([0.75600,0,0])
	coords[0] -= np.array([r-1.512,0,0])
	coords[2:5] -= np.array([r-1.512,0,0])
	atoms = [mol.Atom('C', coords[0]), mol.Atom('C', coords[1])]
	atoms += [mol.Atom('H', c) for c in coords[2::]]
	molecule = mol.Molecule(atoms=atoms)
	molecule.name = f'Ethane {round(r,1)}'
	molecule.save_to_xyz(comment=f'{ff.get_energy(molecule)} {r}')
