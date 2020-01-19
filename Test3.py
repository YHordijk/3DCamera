from modules.screen2 import *
from modules.shape import *
import modules.renderer as rend 
import modules.colour_maps as cmap
import math
import numpy as np
import random
import os

#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = Screen3D(SIZE, camera_position=(4., 4., 10.), camera_orientation=(-0.35,0.35,0), bkgr_colour=(0,0,0))
renderer = rend.Renderer(SIZE, colour_map=cmap.CoolWarm())
clock = pg.time.Clock()
FPS = 100
run = True


#main loop

tick = clock.tick_busy_loop
updt = 0
time = 0


config = [1,1,1,1]


class GTAO:
	def __init__(self, ao_type, n=5):
		self.ao_type = ao_type
		self.coords = 0

		self.n = n
		self.l = list(range(n+1))[::-2]
		self.k = [k//2 for k in list(range(n+1))[::2]]

	@staticmethod
	def get_func(ao_type):
		if ao_type == '1s':
			return lambda x, y, z: (0.08724638 * np.exp(-0.151623*(x**2+y**2+z**2)) + 0.27181242724*np.exp(-0.851819*(x**2+y**2+z**2)))
		elif ao_type == '2s':
			return lambda x, y, z: 0.61282*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2px':
			return lambda x, y, z: 0.61282*x*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*x*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2py':
			return lambda x, y, z: 0.61282*y*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*y*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2pz':
			return lambda x, y, z: 0.61282*z*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*z*np.exp(-0.851819*(x**2+y**2+z**2))
		return 

orb_1s = GTAO.get_func('1s')
orb_2s = GTAO.get_func('2s')
orb_2px = GTAO.get_func('2px')
orb_2py = GTAO.get_func('2py')
orb_2pz = GTAO.get_func('2pz')
maxi = 0 
mapper = lambda x: x/maxi * 255
mapper = cmap.BlackWhite()

array = np.zeros(SIZE)


# x_range = np.linspace(-8,8,5)
# y_range = np.linspace(-8,8,5)
# z_range = np.linspace(-8,8,5)

# x, y, z = np.meshgrid(x_range, y_range, z_range)

# print(x)
# d = (orb_2py(x+3, y, z) + orb_2py(x-3, y, z))**2
# x, y = screen.project_array((x, y, z))

# array[x, y] = d

# print(array)

# screen.blit_array(array)

# screen.draw_axes(2)
# screen.show()

draw_circle = screen.draw_circle

file = os.getcwd() + f'\\Molecules\\aspirin.xyz'
mol = np.loadtxt(file, skiprows=2, usecols=(1,2,3), dtype=float)
print(mol)
samples = 10_000_0
points = 20000
x, y, z = ((np.random.randint(-80000, 80000, size=samples)/10000), (np.random.randint(-80000, 80000, size=samples)/10000), (np.random.randint(-80000, 80000, size=samples)/10000))

d = np.zeros((samples))
for atom in mol:
	d += orb_1s(x-atom[0], y-atom[1], z-atom[2])

# d = orb_1s(x-mol[0][0], y-mol[0][1], z-mol[0][2])

d = d**2




# samples = 10_000_00
# points = 20000

# x, y, z = ((np.random.randint(-80000, 80000, size=samples)/10000), (np.random.randint(-80000, 80000, size=samples)/10000), (np.random.randint(-80000, 80000, size=samples)/10000))
# # d = (orb_2py(x+3, y, z) + orb_2py(x-3, y, z))**2
# d = (config[0]*orb_2s(x,y,z) + config[1]*orb_2py(x,y,z) + config[2]*orb_2px(x,y,z) + config[3]*orb_2pz(x,y,z))**2
# d = orb_2pz(x,y,z)**2
maxi = np.amax(d)
print(d)
colours = mapper[d].T
print(colours)
index = np.arange(0, samples)
index = np.where(abs(d) > maxi/15, index, 0)
index = index[index > 0]

x, y, z, d, colours = x[index][0:points], y[index][0:points], z[index][0:points], d[index][0:points], colours[index][0:points]
# x, y, z, d = x[index][0:points], y[index][0:points], z[index][0:points], d[index][0:points]


# colours = np.array([mapper[d], mapper[d], mapper[d]]).astype(int).T


maxc = math.sqrt(abs(np.amax(x)) + abs(np.amax(y)))*3
print(maxc)

# circles = []
# while len(circles) < 10000:
# 	for _ in range(10000):
# 		x, y, z = (random.randrange(-80000, 80000)/10000,random.randrange(-80000, 80000)/10000,random.randrange(-80000, 80000)/10000)
# 		# d = (orb_2py(x+3, y, z) - orb_2px(x-3, y, z))**2

# 		# d = (orb_2py(x+3, y, z) - orb_2py(x-3, y, z))**2
# 		d = (orb_2s(x,y,z) + orb_2py(x,y,z) + orb_2px(x,y,z))**2
# 		if d > maxi: 
# 			maxi = d
# 			circles = []
# 			break
# 		circles.append((np.asarray((x, y, z)), int(math.sqrt(mapper(d))), (mapper(d),mapper(d),mapper(d))))



while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()

	screen.clear()

	

	screen.camera_position = np.array([maxc*sin(math.pi * time/5), 3, maxc*cos(math.pi * time/5)])
	screen.camera_orientation = np.array([-.3, math.pi * time/5, 0])
	
	#code


	# d = (orb_2py(x+3, y, z) + orb_2py(x-3, y, z))**2

	# x, y = screen.project_array((x, y, z))

	# x, y = d

	# x, y, z = (np.random.randint(-80000, 80000, size=50000)/10000, np.random.randint(-80000, 80000, size=50000)/10000, np.random.randint(-80000, 80000, size=50000)/10000)
	# d = (orb_2py(x+3, y, z) + orb_2py(x-3, y, z))**2


	# pos = np.asarray(np.vstack((x, y, z)).T)
	# selection = np.random.choice(np.arange(50000), size=10000, p=d/sum(d))
	# pos = pos[selection]
	# d = d[selection]

	# pos = screen.project_array(pos)

	# x, y  = np.hsplit(pos, 2)
	# x, y = x.flatten(), y.flatten()
	# x = np.minimum(x, WIDTH-1)
	# y = np.minimum(y, HEIGHT-1)

	# array[x,y] = d

	for atom in mol:
		draw_circle(atom, 5, (255,0,0))

		# screen.draw_pixels(np.asarray((x, y, z)), (mapper(d),mapper(d),mapper(d)))

	



	screen.draw_pixels(np.asarray((x, y, z)).T, colour_array=colours)

	screen.draw_axes(4)
	# tick end
	screen.update()


	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
	

# renderer.input_array(array)
# renderer.show()

