from modules.screen2 import *
from modules.shape import *
import math
import numpy as np
import random

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

orb_1s = GTAO.get_func('1s')
orb_2px = GTAO.get_func('2px')
orb_2py = GTAO.get_func('2py')

def Lorentz(pos, dT, sigma=10, ro=28, beta=8/3):
	x, y, z = pos
	x += dT * (sigma * (pos[1] - pos[0]))
	y += dT * (pos[0] * (ro - pos[2]) - pos[1])
	z += dT * (pos[0] * pos[1] - beta * pos[2])
	
	return np.array([x, y, z])

def Aizawa(pos, dT, a=0.95, b=0.7, c=0.6, d=3.5, e=0.25, f=0.1):
	x, y, z = pos
	nx = x + dT * ((z-b)*x - d*y)
	ny = y + dT * (d*x + (z-b)*y)
	nz = z + dT * (c + a*z - z**3/3 - (x**2+y**2) * (1+e*z) + f*z*x**3)
	
	return np.array([nx, ny, nz])

def Bouali (pos, dT, a=0.3, b=1.):
	x, y, z = pos
	nx = x + dT * (x*(4-y) + a*z)
	ny = y + dT * (-y*(1-x**2))
	nz = z + dT * (-x*(1.5-b*z)-0.05*z)
	
	return np.array([nx, ny, nz])

#game setup
WIDTH, HEIGHT = SIZE = (1600, 900)
s = Screen3D(SIZE, camera_position=(0., 0., 20.))
clock = pg.time.Clock()
FPS = 120
run = True



#lorentz parameters
curr_pos = np.array([1,1,1])

#main loop
updt = 0
time_move = 0

skip_frames = 1

poss = []

tick = clock.tick_busy_loop
cam_pos = s.camera_position
cam_orn = s.camera_orientation

while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000 / 8
	
	keys = pg.key.get_pressed()

	
	if pg.mouse.get_pressed()[2]:
		pg.event.set_grab(True)
		pg.mouse.set_visible(False)
		pg.mouse.set_pos = (WIDTH//2, HEIGHT//2)
		updt = 0
		time_move += 1

		s.clear()
		
		if keys[pg.K_a]:
			cam_pos[0] += dT * 100
		if keys[pg.K_d]:
			cam_pos[0] -= dT * 100
		if keys[pg.K_w]:
			# cam_pos[2] -= dT * 100
			s.move(((-dT * 100), 0, 0))
		if keys[pg.K_s]:
			# cam_pos[2] += dT * 100
			s.move(((dT * 100), 0, 0))
		if keys[pg.K_SPACE]:
			cam_pos[1] += dT * 100
		if keys[pg.K_LCTRL]:
			cam_pos[1] -= dT * 100	

		move = pg.mouse.get_rel()
		if time_move > 5:
			cam_orn[0] += -(move[1])/500
			cam_orn[1] += (move[0])/300
			cam_orn[0] = max(-0.5*math.pi, min(cam_orn[0], 0.5*math.pi))

		if len(poss) > 500:
			poss = poss[-500:]

		s.draw_axes(10)


	elif updt == 1:
		time_move = 0
		pg.event.set_grab(False)
		pg.mouse.set_visible(True)
		s.clear()


	#code
	
	poss.append(curr_pos)
	for i in range(skip_frames):
		curr_pos = Bouali(curr_pos, dT)

	if len(poss) > 1000:
		poss.pop(0)

	try:
		# for p in poss:
		# 	s.draw_pixel(p)
		s.draw_pixels(poss)
	except:
		raise

	#tick end
	s.update()
	if keys[pg.K_ESCAPE]:
		run = False

	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False