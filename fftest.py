import modules.molecule6 as mol 
import modules.reaxff as reaxff
import modules.uff as uff
import modules.utils as utils
import os
import numpy as np

utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(True)



# molecule = mol.Molecule('ethane')
# # coords = np.asarray([a.coords for a in molecule.atoms])

# # coords -= np.array([0.75600,0,0])
# # coords[0] += np.array([0.512,0,0])
# # coords[2:5] += np.array([0.512,0,0])
# # atoms = [mol.Atom('C', coords[0]), mol.Atom('C', coords[1])]
# # atoms += [mol.Atom('H', c) for c in coords[2::]]
# # molecule = mol.Molecule(atoms=atoms)


# print(molecule)


# ff = reaxff.ForceField(molecule)

# print(np.round(ff.bond_orders, 2))

molecule = mol.Molecule('glycine')
ff = uff.ForceField()

print(ff.get_energy(molecule))