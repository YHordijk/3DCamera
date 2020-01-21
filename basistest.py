import modules.basisset6 as bs 
import modules.molecule6 as mol
import modules.screen4 as scr
import numpy as np 
import pygame as pg
import os
from math import sin, pi, cos

molecule = mol.Molecule('propane.pcp')


# atoms = [mol.Atom('H', (0,0,0)), mol.Atom('H', (1,0,0))]
# atoms = [mol.Atom('H', (0, 1.43233673, -0.96104039)), mol.Atom('H', (0, -1.43233673, -0.96104039)), mol.Atom('O', (0, 0, 0.24026010))]
# atoms = [mol.Atom('C', (0, 0, 0)), mol.Atom('H', (0, 0, 1.1)), mol.Atom('H', (1.03709, 0, -0.366667)), mol.Atom('H', (-0.518545, 0.898146, -0.366667)), mol.Atom('H', (-0.518545, -0.898146, -0.366667))]



# molecule = mol.Molecule(atoms=atoms, basis_set_type='STO-6G')

bs.extended_huckel(molecule)

mos = molecule.molecular_orbitals

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


while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()

	screen.clear()

	screen.camera_position = np.array([5*sin(pi * time/5), 2, 5*cos(pi * time/5)])
	screen.camera_orientation = np.array([-.3, pi * time/5, 0])
	

	screen.draw_axes(1)
	screen.draw_density(mos[11], points=20000)
	screen.draw_shape(molecule, wireframe=True)



	
	# tick end
	screen.update()


	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
	

