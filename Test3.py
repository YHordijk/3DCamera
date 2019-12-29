from modules.screen2 import *
from modules.shape import *
import modules.platformfuncs as pf
import math
import numpy as np
import random

#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = Screen3D(SIZE, camera_position=(0., 0., 10.), camera_orientation=(0,0,0))
clock = pg.time.Clock()
FPS = 100
run = True


player = pf.Player(position=(0.,5.,0.), length=5, width=5, height=7, centering="center bottom")
floor = Rectangle(position=(0.,0.,0.), length=10, width=50, height=30, centering="center top")
p_vel = player.velocity
np_norm = np.linalg.norm


#main loop

tick = clock.tick_busy_loop

updt = 0
time = 0


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
			return lambda x, y, z: 0.08724638*(x**3 + y**3 + z**3) * np.exp(-0.151623*(x**2+y**2+z**2)) + 0.27181242724*(x**3 + y**3 + z**3)*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2s':
			return lambda x, y, z: 0.61282*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2px':
			return lambda x, y, z: 0.61282*x*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*x*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2py':
			return lambda x, y, z: 0.61282*y*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*y*np.exp(-0.851819*(x**2+y**2+z**2))
		elif ao_type == '2pz':
			return lambda x, y, z: 0.61282*z*np.exp(-0.151623*(x**2+y**2+z**2)) + 0.0494718*z*np.exp(-0.851819*(x**2+y**2+z**2))
		return 

orb_2s = GTAO.get_func('2s')
orb_2px = GTAO.get_func('2px')
orb_2py = GTAO.get_func('2py')

mapper = lambda x: x/0.9 * 255

maxi = 0 
while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()
	# screen.clear()
	screen.draw_axes(2)


	#code
	# x_range = np.linspace(-8,8,1600)
	# y_range = np.linspace(-8,8,1600)
	# z_range = np.linspace(-8,8,1600)

	# x, y, z = np.meshgrid(x_range, y_range, z_range)

	# d = (orb_2py(x+3, y, z) + orb_2py(x-3, y, z))**2

	# x, y = screen.project_array((x, y, z))

	# x, y = d

	# x, y, z = (np.random.randint(-80000, 80000, size=50000)/10000, np.random.randint(-80000, 80000, size=50000)/10000, np.random.randint(-80000, 80000, size=50000)/10000)
	# d = (orb_2py(x+3, y, z) + orb_2py(x-3, y, z))**2
	# if np.amax(d) > maxi:
	# 	maxi = np.amax(d)
	# 	print(maxi)
	# pos = np.asarray(np.vstack((x, y, z)).T)
	# pos = pos[np.random.choice(np.arange(50000), size=10000, p=d/sum(d))]
	# screen.draw_pixels(pos)

	for _ in range(1000):
		x, y, z = (random.randrange(-80000, 80000)/10000,random.randrange(-80000, 80000)/10000,random.randrange(-80000, 80000)/10000)
		# d = (orb_2py(x+3, y, z) - orb_2px(x-3, y, z))**2
		d = (orb_2px(x+3, y, z) - orb_2px(x-3, y, z))**2
		if d > maxi:
			maxi = d
			print(maxi)

		screen.draw_circle(np.asarray((x, y, z)), int(math.sqrt(mapper(d))), (mapper(d),mapper(d),mapper(d)))

		# screen.draw_pixels(np.asarray((x, y, z)), (mapper(d),mapper(d),mapper(d)))


	#tick end
	screen.update()

	if keys[pg.K_ESCAPE]:
		run = False
		print(player.position)
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break




