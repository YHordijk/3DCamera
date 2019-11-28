import modules.screen3 as scr 
import modules.molecule as mol 
import pygame as pg
import math, os
import numpy as np
from time import perf_counter
import sys
import pubchempy as pcp

input_mol = 'Aminoindane'



# if len(sys.argv) == 1:
# 	m = mol.Molecule(file=os.getcwd() + rf'\Molecules\{input_mol}.xyz')
# else:
# 	print(os.path.exists(sys.argv[1]))
# 	m = mol.Molecule(file=sys.argv[1])

m = mol.Molecule(molecule_file='Beta-carotene.pcp')

# m.remove_by_element('H')
# m.add_hydrogens(1.1)
m.center()

#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=(100, 100, 190))
clock = pg.time.Clock()
FPS = 120
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
	screen.draw_shape(m, draw_atoms=True)

	keys = pg.key.get_pressed()

	if keys[pg.K_SPACE]:
		m.position = [0,0,0]

	move = pg.mouse.get_rel()
	if pg.mouse.get_pressed()[2]:
		m.position[0] -= move[0]/50
		m.position[1] -= move[1]/50

	if pg.mouse.get_pressed()[0]:
		if keys[pg.K_LCTRL] or keys[pg.K_RCTRL]:
			rot = np.asarray([0, 0, (abs(move[0]) + abs(move[1]))/250])
		else:
			rot = np.asarray([move[1]/250, -move[0]/150, max(-0.5*math.pi, min(rot[0], 0.5*math.pi))])
	else:
		rot *= 0.8
	m.rotate(rot)

	event = pg.event.get(eventtype=pg.MOUSEBUTTONDOWN)

	if len(event) > 0:
		if event[0].button == 4:
			if screen.camera_position[2] > 5:
				screen.camera_position[2] -= dT * screen.camera_position[2] * 3
		if event[0].button == 5:
			if screen.camera_position[2] < 40:
				screen.camera_position[2] += dT * screen.camera_position[2] * 3

	#tick end
	screen.update()

	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
