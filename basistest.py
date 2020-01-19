import modules.basisset2 as bs 
import modules.molecule4 as mol
import numpy as np 
import os

# molecule = mol.Molecule(os.getcwd() + f'\\Molecules\\altrose.xyz', basis_set_type='STO-6G')
# molecule = mol.Molecule('fluorine.pcp', basis_set_type='STO-6G')

# print(molecule.atoms[0].atomic_orbitals.primitives_list)
# print(molecule.atoms)


g1 = bs.AtomicOrbital('STO-6G', mol.Atom('C', (0,0,0)))
g2 = bs.AtomicOrbital('STO-6G', mol.Atom('H', (0,0,1)))

print(g1.primitives_list)
print(g2.primitives_list)

print(bs.overlap_integral(g1.primitives_list[2], g2.primitives_list[0]))




# print(bs.extended_huckel(molecule))

