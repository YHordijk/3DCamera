import modules.screen2 as scr 
import modules.shape as shp 
import pygame as pg
import math, os
import numpy as np

m = shp.Molecule(position=(0,0,0))
m.load_xyz(os.getcwd() + f'\Molecules\hexabenzocoronene2.xyz')
# m.load_xyz(r'C:\Users\Yuman\Desktop\Programmeren\Python\school\Introduction_to_scientific_programming\Week_5\molecules\benzene_with_hydrogens.xyz')


# m2 = shp.Molecule(position=(-5,0,0))
# m2.load_xyz(os.getcwd() + r'\Molecules\pyridine.xyz')


#game setup
WIDTH, HEIGHT = SIZE = (1600, 900)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0))
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
	
	screen.clear()

	screen.draw_shape(m)

	# d = np.sqrt(sum((screen.camera_position - m.position)**2))
	# d2 = np.sqrt(sum((screen.camera_position - m2.position)**2))

	# if d >= d2:
	# 	screen.draw_shape(m)
	# 	screen.draw_shape(m2)
	# else:
	# 	screen.draw_shape(m2)
	# 	screen.draw_shape(m)


	# m.rotation = (0, math.sin(time/10)/50, math.sin(time/10)/50)
	# m.position = (5*math.sin(time), 0, 5*math.cos(time))
	# m2.rotation = (0, -math.sin(time/10+1.6)/50, -math.sin(time/10+1.6)/50)
	# m2.position = (-5*math.sin(time), 0, -5*math.cos(time))
	
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
