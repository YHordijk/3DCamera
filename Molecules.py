import modules.screen3 as scr 
import modules.molecule as mol 
import pygame as pg
import math, os
import numpy as np
from time import perf_counter
import sys


if sys.argv is None:
	m = mol.Molecule(file=os.getcwd() + r'\Molecules\Fullerene.xyz')
else:
	print(os.path.exists(sys.argv[1]))
	m = mol.Molecule(file=sys.argv[1])
# m.remove_hydrogens()
# m.add_hydrogens(1.1)
m.center()


#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=(0,0,0))
clock = pg.time.Clock()
FPS = 60
run = True

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0

rot = np.array([0.,0.,0.])

while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	
	screen.clear()
	screen.draw_shape(m)

	move = pg.mouse.get_rel()
	if pg.mouse.get_pressed()[0]:
		rot = np.asarray([move[1]/250, -move[0]/150, max(-0.5*math.pi, min(rot[0], 0.5*math.pi))])
	else:
		rot *= 0.6
	m.rotate(rot)

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
