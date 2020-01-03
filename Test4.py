from modules.screen2 import *
from modules.shape import *
import modules.molecule4 as mol
import modules.renderer as rend 
import modules.colour_maps as cmap
import math
import numpy as np
import random
import os

#game setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = Screen3D(SIZE, camera_position=(4., 4., 10.), camera_orientation=(-0.35,0.35,0), bkgr_colour=(0,0,0))
renderer = rend.Renderer(SIZE, colour_map=cmap.CoolWarm())
clock = pg.time.Clock()
FPS = 100
run = True


#main loop

tick = clock.tick_busy_loop
updt = 0
time = 0

molecule = mol.Molecule(os.getcwd() + f'\\Molecules\\ethane.xyz', basis_set_type='STO-6G')

rang = 4
samples = 1000000
x, y, z = ((np.random.randint(-rang*10000, rang*10000, size=samples)/10000), (np.random.randint(-rang*10000, rang*10000, size=samples)/10000), (np.random.randint(-rang*10000, rang*10000, size=samples)/10000))
d = molecule.get_orb_density(np.asarray((x, y, z)).T).flatten()

index = np.arange(0, samples)
print(np.amax(d)/15)
index = np.where(abs(d) > 20, index, 0)
index = index[index > 0]

mapper = cmap.BlackWhite()
colours = mapper[d].T

points = 50000
x, y, z, d, colours = x[index][0:points], y[index][0:points], z[index][0:points], d[index][0:points], colours[index][0:points]


maxc = 4

draw_circle = screen.draw_circle

while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()

	screen.clear()

	

	screen.camera_position = np.array([maxc*sin(math.pi * time/5), 1.5, maxc*cos(math.pi * time/5)])
	screen.camera_orientation = np.array([-.3, math.pi * time/5, 0])
	

	# for atom in molecule.atoms:
	# 	draw_circle(atom.coords, 5, (255,0,0))


	screen.draw_pixels(np.asarray((x, y, z)).T, colour_array=colours)

	screen.draw_axes(0.2)
	# tick end
	screen.update()


	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break
	

# renderer.input_array(array)
# renderer.show()

