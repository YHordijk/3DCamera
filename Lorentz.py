from modules.screen2 import *
from modules.shape import *
import math
import numpy as np
import random


#game setup
WIDTH, HEIGHT = SIZE = (600, 600)
s = Screen3D(SIZE, camera_position=(0, 0, 120))
clock = pg.time.Clock()
FPS = 120
run = True

pg.event.set_grab(True)
pg.mouse.set_visible(False)

#lorentz parameters
curr_pos = np.array([1,1,1])
sigma = 10
ro = 28
beta = 8/3

#main loop
updt = 0



poss = []

while run:
	#tick prep
	updt += 1
	dT = clock.tick_busy_loop(FPS)/1000
	s.clear()

	pg.mouse.set_pos = (500, 300)

	keys = pg.key.get_pressed()
	if keys[pg.K_LEFT]:
		s.camera_position[0] += 1
	if keys[pg.K_RIGHT]:
		s.camera_position[0] -= 1
	if keys[pg.K_UP]:
		s.camera_position[2] -= 1
	if keys[pg.K_DOWN]:
		s.camera_position[2] += 1
	if keys[pg.K_SPACE]:
		s.camera_position[1] += 1
	if keys[pg.K_LCTRL]:
		s.camera_position[1] -= 1

	

	move = pg.mouse.get_rel()
	s.camera_orientation[0] += -(move[1])/500
	s.camera_orientation[1] += (move[0])/300
	s.camera_orientation[0] = max(-0.5*math.pi, min(s.camera_orientation[0], 0.5*math.pi))
	#code
	poss.append(curr_pos)
	s.draw_axes(10)
	x, y, z = curr_pos

	x += dT * (sigma * (curr_pos[1] - curr_pos[0]))
	y += dT * (curr_pos[0] * (ro - curr_pos[2]) - curr_pos[1])
	z += dT * (curr_pos[0] * curr_pos[1] - beta * curr_pos[2])
	
	curr_pos = np.array([x, y, z])

	if len(poss) > 400:
		poss.pop(0)

	try:
		s.draw_lines(poss, closed=False)
	except:
		pass

	#tick end
	s.update()
	if keys[pg.K_ESCAPE]:
		run = False

	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False