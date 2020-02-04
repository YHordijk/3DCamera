import modules.molecule6 as mol 
import modules.reaxff as reaxff
import modules.uff as uff
import modules.utils as utils
import modules.plot2 as plot
import os
import numpy as np

utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(True)

p = plot.Plot()

ff = uff.ForceField()
ethane = mol.Molecule('ethane')
ethane_coords = np.asarray([a.coords for a in ethane.atoms])

x = []
y = []
for r in np.arange(1, 2, 0.1):
	coords = ethane_coords.copy()
	coords -= np.array([0.75600,0,0])
	coords[0] -= np.array([r-1.512,0,0])
	coords[2:5] -= np.array([r-1.512,0,0])
	atoms = [mol.Atom('C', coords[0]), mol.Atom('C', coords[1])]
	atoms += [mol.Atom('H', c) for c in coords[2::]]
	molecule = mol.Molecule(atoms=atoms)

	x.append(r)
	y.append(ff.get_energy(molecule))


p.plot(x,y, style='scatter')
p.y_label = 'Energy (kJ/mol)'
p.x_label = 'r (angstrom)'
p.title = 'Energy as function of length of H3C-CH3 bond'
p.show()