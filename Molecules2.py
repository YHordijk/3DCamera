import modules.screen2 as scr 
import modules.shape as shp 
import pygame as pg
import math, os
from time import perf_counter
import numpy as np

m = shp.Molecule(position=(0,0,0))
m.load_xyz(os.getcwd() + f'\Molecules\hexabenzocoronene2.xyz')
# m.load_xyz(r'C:\Users\Yuman\Desktop\Programmeren\Python\school\Introduction_to_scientific_programming\Week_5\molecules\benzene.xyz')


# m2 = shp.Molecule(position=(-5,0,0))
# m2.load_xyz(os.getcwd() + r'\Molecules\pyridine.xyz')


#game setup
WIDTH, HEIGHT = SIZE = (1600, 900)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=(0,0,0))
clock = pg.time.Clock()
FPS = 30
run = True

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0


while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT

	screen.clear()
	screen.draw_shape(m)

	move = pg.mouse.get_rel()
	if pg.mouse.get_pressed()[0]:
		rot = m.rotation
		m.rotation = [move[1]/250, -move[0]/150, max(-0.5*math.pi, min(rot[0], 0.5*math.pi))]
	else:
		m.rotation = [0.9 *  r for r in m.rotation]

	keys = pg.key.get_pressed()
	if keys[pg.K_UP]:
		if screen.camera_position[2] > 5:
			screen.camera_position[2] -= dT * screen.camera_position[2]
	if keys[pg.K_DOWN]:
		if screen.camera_position[2] < 40:
			screen.camera_position[2] += dT * screen.camera_position[2]

	#tick end
	screen.update()

	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
