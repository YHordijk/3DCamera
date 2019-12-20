import modules.screen3 as scr 
import modules.molecule3 as mol 
import modules.forcefields as ff
import math, os
import numpy as np
from time import perf_counter
import sys
import pubchempy as pcp
import pygame as pg



# mols = [mol.Molecule(molecule_file='ethane.pcp', warning_level=1, scale=400)]

mols = [mol.Molecule(molecule_file=os.getcwd() + f'\\Molecules\\ethane.xyz', warning_level=1, scale=400)]

# mols = [mol.Molecule(molecule_file='Glucose.pcp', warning_level=1, position=[5,0,0], scale=400),
# 		  mol.Molecule(molecule_file='Altrose.pcp', warning_level=1, position=[-5,0,0], scale=400),]


mol = mols[0]
atoms = mol.atoms
# print(mol.torsion_angle(atoms[2], atoms[0], atoms[1], atoms[7]))

# print(id(atoms[0]))

mm = ff.MM2()
print(mm.get_energy(mol))

#screen setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=(100, 100, 190))
pg.display.set_caption(', '.join([m.name.capitalize() for m in mols]))
clock = pg.time.Clock()
FPS = 120
run = True

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0

rot = np.array([0.,0.,0.])
zoom = 0

pg.key.set_repeat()

while run:
	# start = perf_counter()
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT

	screen.draw_axes(100)
	
	screen.clear()
	[screen.draw_shape(m, draw_atoms=True, draw_bonds=True, draw_hydrogens=True) for m in mols]
	[m.rotate(rot) for m in mols]
	# screen.draw_shape(m, draw_atoms=True, draw_bonds=True, draw_hydrogens=True)

	keys = pg.key.get_pressed()

	# if keys[pg.K_SPACE]:
	# 	m.position = [0,0,0]

	move = pg.mouse.get_rel()
	if pg.mouse.get_pressed()[2]:
			screen.camera_position[0] += move[0]/50
			screen.camera_position[1] += move[1]/50

	if pg.mouse.get_pressed()[0]:
		if keys[pg.K_LCTRL] or keys[pg.K_RCTRL]:
			rot = np.asarray([0, 0, (abs(move[0]) + abs(move[1]))/250])
		else:
			rot = np.asarray([move[1]/150, -move[0]/150, 0])
	else:
		rot *= 0.8


	# m.rotate(rot)

	event = pg.event.get(eventtype=pg.MOUSEBUTTONDOWN)

	if len(event) > 0:
		if event[0].button == 4:
			if screen.camera_position[2] > 3:
				zoom = -dT * screen.camera_position[2] * 3
		if event[0].button == 5:
			if screen.camera_position[2] < 40:
				zoom = dT * screen.camera_position[2] * 3

	screen.camera_position[2] += zoom
	zoom *= 0.8

	#tick end
	screen.update()

	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break

	# print(perf_counter()-start)