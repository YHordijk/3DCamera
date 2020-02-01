import modules.molecule6 as mol 
import modules.reaxff as reaxff
import modules.utils as utils
import os
import numpy as np

utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(True)

molecule = mol.Molecule('ethane')


coords = np.asarray([a.coords for a in molecule.atoms])
coords -= 0.75600
coords[0] += 0.512
coords[2:5] += 0.512
atoms = [mol.Atom('C', coords[0]), mol.Atom('C', coords[1])]
atoms += [mol.Atom('H', c) for c in coords[2::]]
molecule = mol.Molecule(atoms=atoms)

print(molecule)


ff = reaxff.ForceField(molecule)

print(np.round(ff.bond_orders, 2))


	