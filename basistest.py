import modules.basisset6 as bs 
import modules.molecule6 as mol
import modules.screen4 as scr
import modules.utils as utils
import numpy as np 
import pygame as pg
import os
from math import sin, pi, cos

utils.ff_print_source(False)
utils.ff_use_colours(False)
utils.ff_print_time(True)

molecule = mol.Molecule(os.getcwd() + f'\\Molecules\\water.xyz')


# atoms = [mol.Atom('H', (0,0,0)), mol.Atom('H', (1,0,0))]
# atoms = [mol.Atom('H', (0, 1.43233673, -0.96104039)), mol.Atom('H', (0, -1.43233673, -0.96104039)), mol.Atom('O', (0, 0, 0.24026010))]
# atoms = [mol.Atom('C', (0, 0, 0)), mol.Atom('H', (0, 0, 1.1)), mol.Atom('H', (1.03709, 0, -0.366667)), mol.Atom('H', (-0.518545, 0.898146, -0.366667)), mol.Atom('H', (-0.518545, -0.898146, -0.366667))]



# molecule = mol.Molecule(atoms=atoms, basis_set_type='STO-6G')

bs.extended_huckel(molecule)
mos = molecule.molecular_orbitals

camera_range = max(max(mos[-1].ranges[0], mos[-1].ranges[1]), max(mos[-1].ranges[2], mos[-1].ranges[3])) + 1

#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = scr.Screen3D(SIZE, camera_position=(3., 3., 10.), camera_orientation=(-0.35,0.35,0), bkgr_colour=(0,0,0))
clock = pg.time.Clock()
FPS = 100
run = True

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0

mo_numb = 0

while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()

	screen.clear()

	screen.camera_position = np.array([camera_range*sin(pi * time/5), camera_range/3, camera_range*cos(pi * time/5)])
	# screen.camera_orientation = np.array([-.3, pi * time/5, 0])
	screen.look_at((0,0,0))

	screen.draw_axes(1)
	screen.draw_density(mos[mo_numb%len(mos)], points=5000)
	screen.draw_shape(molecule, wireframe=True)


	
	# tick end
	screen.update()


	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.KEYDOWN:
			if event.key == pg.K_RIGHT:
				mo_numb += 1
			if event.key == pg.K_LEFT:
				mo_numb -= 1

		if event.type == pg.QUIT:
			run = False
			break
	

