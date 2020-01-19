import modules.basisset3 as bs 
import modules.molecule5 as mol
import modules.renderer as rend
import gbasis.evals.eval as gb_eval
from gbasis.integrals.overlap import overlap_integral
import numpy as np 
import os

# molecule = mol.Molecule('methane.pcp')


# atoms = [mol.Atom('H', (0,0,0)), mol.Atom('H', (1,0,0))]
atoms = [mol.Atom('C', (0,0,0)), mol.Atom('H', (1,0,0)), mol.Atom('H', (1,0,0))]

basis = bs.BasisSet(atoms, 'STO-6G')
print(basis.basis)

# print(len(molecule.basis.basis))

size = 600
grid = np.linspace(-2,2,size)
X, Y = np.meshgrid(grid, grid)
grid_3d = np.vstack([X.ravel(), Y.ravel(), np.zeros(size*size).ravel()]).T

transform = [[0,0,1,1]]
transform = np.asarray([np.asarray(t/np.linalg.norm(t)**2) for t in transform])


a = gb_eval.evaluate_basis(basis.basis, grid_3d, coord_type='cartesian')[2]
a = np.reshape(a, (size,size))

# print(overlap_integral(basis.basis, coord_type='cartesian', transform=transform))


renderer = rend.Renderer()
renderer.input_array(a)
renderer.show()