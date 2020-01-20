import modules.basisset2 as bs 
import modules.molecule5 as mol
import modules.renderer as rend
import numpy as np 
import os

# molecule = mol.Molecule('methane.pcp')


# atoms = [mol.Atom('H', (0,0,0)), mol.Atom('H', (1,0,0))]
atoms = [mol.Atom('H', (0, 1.43233673, -0.96104039)), mol.Atom('H', (0, -1.43233673, -0.96104039)), mol.Atom('O', (0, 0, 0.24026010))]



molecule = mol.Molecule(atoms=atoms, basis_set_type='STO-3G')

# basis = bs.BasisSet('STO-3G', atoms)
bs.extended_huckel(molecule)



# basis.basis.get_integrals()
# print(basis.basis.overlap_integrals)

# # print(len(molecule.basis.basis))

# size = 600
# grid = np.linspace(-2,2,size)
# X, Y = np.meshgrid(grid, grid)
# grid_3d = np.vstack([X.ravel(), Y.ravel(), np.zeros(size*size).ravel()]).T

# transform = [[0,0,1,1]]
# transform = np.asarray([np.asarray(t/np.linalg.norm(t)**2) for t in transform])


# a = gb_eval.evaluate_basis(basis.basis, grid_3d, coord_type='cartesian')[2]
# a = np.reshape(a, (size,size))

# # print(overlap_integral(basis.basis, coord_type='cartesian', transform=transform))


# renderer = rend.Renderer()
# renderer.input_array(a)
# renderer.show()