import math
import numpy as np
import modules.renderer as rend
import modules.colour_maps as cmap
import scipy.integrate

class GTAO:
	def __init__(self, ao_type='2px', n=5):
		self.ao_type = ao_type
		self.coords = 0

		self.n = n
		self.l = list(range(n+1))[::-2]
		self.k = [k//2 for k in list(range(n+1))[::2]]

	@staticmethod
	def get_func(ao_type):
		if ao_type == '1s':
			return lambda x, y, z: (0.08724638* np.exp(-0.151623*(x**2+y**2+z**2)) + 0.27181242724*np.exp(-0.851819*(x**2+y**2+z**2)))/math.sqrt(10.15350166792534)
		elif ao_type == '2s':
			return lambda x, y, z: (0.61282*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*np.exp(-0.851819*(x**2+y**2+z**2)))
		elif ao_type == '2px':
			return lambda x, y, z: (0.61282*x*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*x*np.exp(-0.851819*(x**2+y**2+z**2)))/math.sqrt(10.49388878806781*2)
		elif ao_type == '2py':
			return lambda x, y, z: (0.61282*y*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*y*np.exp(-0.851819*(x**2+y**2+z**2)))/math.sqrt(10.49388878806781*2)
		elif ao_type == '2pz':
			return lambda x, y, z: (0.61282*z*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*z*np.exp(-0.851819*(x**2+y**2+z**2)))/math.sqrt(10.49388878806781*2)
		return 

orb_1s = GTAO.get_func('1s')
orb_2px = GTAO.get_func('2px')
orb_2py = GTAO.get_func('2py')

x_range = np.linspace(-8,8,600)
y_range = np.linspace(-8,8,600)
Z = 0
X, Y = np.meshgrid(x_range, y_range)

renderer = rend.Renderer((400,400), colour_map=cmap.CoolWarm())

# for d in np.linspace(7, 3, 50):
d = 4
# dens = (-orb_2px(X - d/2, Y, Z) + orb_2px(X + d/2, Y, Z) - orb_2py(X - d/2, Y, Z) + orb_2py(X + d/2, Y, Z))
# dens = (orb_2py(X - d/2, Y, Z) + orb_2py(X + d/2, Y, Z))
# dens = (orb_1s(X - d/2, Y, Z) - orb_1s(X + d/2, Y, Z))
# dens = orb_1s(X,Y,Z)**2
r = 0

def orb(Z, Y, X, r):
	return orb_1s(X, Y-r, Z) * orb_1s(X, Y+r, Z)

dens = orb(Z, Y, X, 5)
print(dens.shape)
renderer.input_array(dens)
renderer.show()


print(scipy.integrate.tplquad(orb, -8,8,-8,8,-8,8, args=[(0)]))
# print(scipy.integrate.tplquad(orb_1s, -8,8,-8,8,-8,8))
# dens = orb(Z, Y, X, r)
# # dens[0,0] = 1.644058873447782e-112
# renderer.input_array(dens)
# print(np.amax(dens))
# renderer.show()