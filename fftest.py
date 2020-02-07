import modules.molecule6 as mol 
import modules.reaxff as reaxff
import modules.uff as uff
import modules.utils as utils
import modules.plot as plot
import modules.minimizer as mini
import os, math
import numpy as np

utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(False)



ff = uff.ForceField()

ethane = mol.Molecule('1,2-dichloroethane')
# ethane.rotate_bond(ethane.atoms[10], ethane.atoms[5], dT*math.pi)

print(ethane.atoms)

p = plot.Plot()



# x = []
# y = []
# for r in np.arange(1,3,0.05):
# 	ethane.stretch_bond(ethane.atoms[0], ethane.atoms[1], r)
# 	x.append(r)
# 	y.append(ff.get_energy(ethane, morse_potential=True, verbosity=0))

# p.clear()
# p.plot(x,y, style='line')
# p.y_label = 'Energy (kJ/mol)'
# p.x_label = 'r (angstrom)'
# p.title = 'Energy as function of length of H3C-CH3 bond'
# p.show()

x = []
y = []
ran = 120
for i in range(ran+1):
	ethane.rotate_bond(ethane.atoms[2], ethane.atoms[3], 2*math.pi/ran)
	x.append(i/ran*360)
	y.append(ff.get_energy(ethane, morse_potential=True, verbosity=2))

print(np.argmin(y)*math.pi/180)
p.clear()
p.plot(x,y, style='line')
p.y_label = 'Energy (kJ/mol)'
p.x_label = 'theta (degrees)'
p.title = 'Energy as function of rotation of H3C-CH3 bond'
p.show()


