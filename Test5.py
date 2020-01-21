from modules.screen3 import *
import modules.basisset as bs
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

basis = bs.BasisSet('STO-6G', [mol.Atom('Cl', (0,0,0))])

colour_map = cmap.CoolWarm()

def dens():
	points = 50000
	samples = 50*points
	rang = 4

	x, y, z = ((np.random.randint(-rang*300, rang*300, size=samples)/300), (np.random.randint(-rang*300, rang*300, size=samples)/300), (np.random.randint(-rang*300, rang*300, size=samples)/300))
	d = basis(np.asarray((x, y, z)).T).flatten()**2
	index = np.arange(0, samples)
	index = np.where(abs(d) > np.amax(abs(d))/15, index, 0)

	index = index[index > 0]

	colours = colour_map[d].T

	x, y, z, colours = x[index][0:points], y[index][0:points], z[index][0:points], colours[index][0:points]
	pos = np.asarray((x, y, z)).T
	return pos, colours

pos, colours = dens()


while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()

	screen.clear()

	screen.camera_position = np.array([1*sin(math.pi * time/5), 0.5, 1*cos(math.pi * time/5)])
	screen.camera_orientation = np.array([-.3, math.pi * time/5, 0])
	

	screen.draw_axes(0.2)
	print(colours)
	screen.draw_pixels(pos, colour_array=colours)


	
	# tick end
	screen.update()


	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
	

