import modules.screen4 as scr 
import modules.marchingcube as mc
import numpy as np
import pygame as pg
from math import sin, cos, pi


pg.init()
screen = scr.Screen3D((1200,720))


# array = np.random.random((3,3,3))
# mesh = np.asarray(mc.marching_cubes(array, 0.5, 1)) - 1



array = np.array([[[0,0,0],[0,0.5,0],[0,0,0]],
				  [[0,0.5,0],[0,1,0],[0,0.5,0]],
				  [[0,0,0],[0,0.5,0],[0,0,0]]])

mesh = mc.Mesh(array, 0.5, 1, -1)
# mesh = np.asarray(mc.marching_cubes(array, 0.5, 1)) -1

x = y = z = np.array([-1,0,1])
x, y, z = np.meshgrid(x, y, z)
x, y, z = x.flatten(), y.flatten(), z.flatten()



# rang = 6
# points = 6
# x, y, z = np.linspace(-rang, rang, points), np.linspace(-rang, rang, points), np.linspace(-rang, rang, points)
# x, y, z = np.meshgrid(x, y, z)
# x, y, z = x.flatten(), y.flatten(), z.flatten()

# func = lambda x, y, z: (y**2 - 0.5*z**2 -0.5*x**2)

# array = func(x, y, z).reshape((points,points,points))
# print(array.max())

# mesh = np.asarray(mc.marching_cubes(array, 1, 1)) - rang/2

# print(len(mesh))




clock = pg.time.Clock()
FPS = 120
tick = clock.tick_busy_loop
updt = 0
time = 0
run = True

camera_range = 2


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



	# light = (10*sin(pi * time/5), 10*cos(pi * time/5), 0)
	# light = screen.camera_position
	light = (0,0,10)

	# screen.camera_position = np.array([camera_range*sin(pi * time/5), camera_range/3, camera_range*cos(pi * time/5)])
	# screen.camera_orientation = np.array([-.3, pi * time/5, 0])

	screen.camera_position = np.array([camera_range*sin(pi * 2/5), camera_range/3, camera_range*cos(pi * 2/5)])
	screen.camera_orientation = np.array([-.3, pi * 2/5, 0])



	# screen.draw_pixels(np.asarray((x,y,z)).T)
	screen.draw_axes(1)
	mesh.rotate((0.003,0.,0.))
	screen.draw_mesh(mesh.mesh, lighting=light, fill=False, lighting_colour=(50,200,200))

	# screen.draw_line(np.asarray(((0,0,0), light)), (253, 184, 19))







	screen.update()