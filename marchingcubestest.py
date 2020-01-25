import modules.screen4 as scr 
import modules.marchingcube as mc
import numpy as np
import pygame as pg
from math import sin, cos, pi


pg.init()
screen = scr.Screen3D((1200,720))





# array = np.array([[[0,0,0],[0,0,0],[1,0,1]],
# 				  [[0,2,0],[0,0.75,0],[0,2,0]],
# 				  [[0,0,0],[0,0,0],[1,0,1]]])

# x = y = z = np.array([-1,0,1])
# x, y, z = np.meshgrid(x, y, z)
# x, y, z = x.flatten(), y.flatten(), z.flatten()


# triangles = np.asarray(mc.marching_cubes(array, 0.5, 1)) -1


rang = 6
points = 6
x, y, z = np.linspace(-rang, rang, points), np.linspace(-rang, rang, points), np.linspace(-rang, rang, points)
x, y, z = np.meshgrid(x, y, z)
x, y, z = x.flatten(), y.flatten(), z.flatten()

func = lambda x, y, z: (y**2 - 0.5*z**2 -0.5*x**2)

array = func(x, y, z).reshape((points,points,points))
print(array.max())

triangles = np.asarray(mc.marching_cubes(array, 1, 1)) - rang/2

print(len(triangles))




clock = pg.time.Clock()
FPS = 120
tick = clock.tick_busy_loop
updt = 0
time = 0
run = True

camera_range = rang+1


while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT
	keys = pg.key.get_pressed()
	ev = pg.event.get()
	

	if keys[pg.K_ESCAPE]:
		run = False

	for e in ev:
		if e.type == pg.QUIT:
			run = False

	screen.clear()





	screen.camera_position = np.array([camera_range*sin(pi * time/5), camera_range/3, camera_range*cos(pi * time/5)])
	screen.camera_orientation = np.array([-.3, pi * time/5, 0])
	# screen.draw_pixels(np.asarray((x,y,z)).T)
	screen.draw_mesh(triangles)







	screen.update()