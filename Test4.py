from modules.screen3 import *
from modules.shape import *
import modules.molecule4 as mol
import modules.colour_maps as cmap
import math
import numpy as np
import random
import os

#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = Screen3D(SIZE, camera_position=(4., 4., 10.), camera_orientation=(-0.35,0.35,0), bkgr_colour=(0,0,0))
clock = pg.time.Clock()
FPS = 100
run = True

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0

molecule = mol.Molecule(os.getcwd() + f'\\Molecules\\methane.xyz', basis_set_type='STO-6G')

while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()

	screen.clear()

	screen.camera_position = np.array([4*sin(math.pi * time/5), 1.5, 4*cos(math.pi * time/5)])
	screen.camera_orientation = np.array([-.3, math.pi * time/5, 0])
	

	screen.draw_axes(1)
	screen.draw_density(molecule, 15000, colour_map=cmap.CoolWarm(), mo=6)
	screen.draw_shape(molecule, wireframe=True)

	
	# tick end
	screen.update()


	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
	
