import modules.screen2 as screen
import modules.shape as shape
from math import cos, sin
import math
import numpy as np
import random
import pygame as pg




s = screen.Screen3D((1000,600), camera_position=(10,0,100))
cube = shape.Cube(20, position=(0, 0, 0))
clock = pg.time.Clock()

pg.event.set_grab(True)
pg.mouse.set_visible(False)

final = True
FPS = 30
time = 0
c = 0

while final:
	dT = clock.tick_busy_loop(FPS)/1000
	time += dT
	c += 1
	s.clear()

	pg.mouse.set_pos = (500, 300)

	keys = pg.key.get_pressed()
	if keys[pg.K_LEFT]:
		s.camera_position[0] += dT*50
	if keys[pg.K_RIGHT]:
		s.camera_position[0] -= dT*50
	if keys[pg.K_UP]:
		s.camera_position[2] -= dT*50
	if keys[pg.K_DOWN]:
		s.camera_position[2] += dT*50
	if keys[pg.K_SPACE]:
		s.camera_position[1] += dT*50
	if keys[pg.K_LCTRL]:
		s.camera_position[1] -= dT*50

	

	move = pg.mouse.get_rel()
	s.camera_orientation[0] += -(move[1])/500
	s.camera_orientation[1] += (move[0])/300
	s.camera_orientation[0] = max(-0.5*math.pi, min(s.camera_orientation[0], 0.5*math.pi))

	s.draw_axes(5)
	m = 125
	s.draw_shape(cube, colour=(m + sin((time)*9)*m, 255, 255 - (m + sin((time)*9)*m)))

	cube.rotation[1] += 2*dT
	cube.rotation[0] += dT
	cube.rotation[2] += dT/2

	for j in range(-20, 20):
		s.draw_lines([3*np.array([j*cos(i/2.5), j, j*sin(i/2.5)]) for i in range(0,16)], colour=(m + sin((time + j*math.pi/5)*9)*m, 255, 255 - (m + sin((time + j*math.pi/5)*9)*m)))

	s.update()

	if keys[pg.K_ESCAPE]:
		final = False

	for event in pg.event.get():
		if event.type == pg.QUIT:
			final = False