import modules.molecule6 as mol 
import modules.reaxff as reaxff
import modules.uff as uff
import modules.utils as utils
import modules.plot as plot
import modules.minimizer as mini
import os, math
import numpy as np

from subprocess import check_output


def printer(m):
	for a1, a2 in m.get_unique_bonds():
		print(a1, a2, a1.distance_to(a2))

	for a1, a2, a3, angle in m.get_unique_bond_angles():
		print(a1,a2,a3,angle)


utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(False)


ff = uff.ForceField()
ethane = mol.Molecule('water')
a = ethane.atoms
# ethane2 = mol.Molecule('benzene2')

p = plot.Plot()


# e1 = ff.get_energy(ethane,verbosity=1)*4.18
# e2 = ff.get_energy(ethane2,verbosity=1)*4.18
# print((e2-e1)/0.001)
# print((44.67307-44.64240)/0.001)
ethane.stretch_bond(a[0], a[1], 3)
newmol = mini.minimize(ethane, ff, steps=1500)
printer(ethane)
# newmol.center()
printer(newmol)

# x = []
# y = []
# y_ob = []
# ran = 120
# for i in range(ran+1):
# 	ethane.rotate_bond(ethane.atoms[0], ethane.atoms[1], 2*math.pi/ran)
# 	x.append(i/ran*360)
# 	y.append(ff.get_energy(ethane, morse_potential=True, verbosity=0))

# 	file = f'C:\\Users\\Yuman\\Desktop\\Programmeren\\Python\\PyGame\\3DCamera\\Molecules\\ethane{int(x[-1])}.xyz'
# 	ethane.save_to_xyz(file=file, comment=f'{x[-1]}')

# 	out = check_output(f'obenergy -ff UFF {file}', shell=True)
# 	y_ob.append(float(out[out.find(b'TOTAL ENERGY')+15:-9]))


# x = []
# y = []
# y_ob = []
# ran = 120
# for r in np.arange(0.5,3,0.05):
# 	ethane.stretch_bond(ethane.atoms[0], ethane.atoms[1], r)
# 	print(r)
# 	x.append(r)
# 	y.append(ff.get_energy(ethane, morse_potential=True, verbosity=1))

# 	file = f'C:\\Users\\Yuman\\Desktop\\Programmeren\\Python\\PyGame\\3DCamera\\Molecules\\ethane{round(r,2)}.xyz'
# 	ethane.save_to_xyz(file=file, comment=f'{round(r,1)}')

# 	out = check_output(f'obenergy -ff UFF {file}', shell=True)
# 	y_ob.append(float(out[out.find(b'TOTAL ENERGY')+15:-9]))





# p.clear()

# p.plot(x,y, style='line')
# p.plot(x,y_ob, style='line')
# p.plot(x,np.asarray(y_ob) - np.asarray(y), style='line')
# p.y_label = 'Energy (kJ/mol)'
# p.x_label = 'theta (degrees)'
# p.title = 'Energy as function of rotation of H3C-CH3 bond'
# p.show()


