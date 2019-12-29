from modules.screen2 import *
from modules.shape import *
import modules.platformfuncs as pf
import math
import numpy as np
import random

#game setup
WIDTH, HEIGHT = SIZE = (1600, 900)
screen = Screen3D(SIZE, camera_position=(0., 6., 20.), camera_orientation=(-.3,0,0))
clock = pg.time.Clock()
FPS = 100
run = True


player = pf.Player(position=(0.,5.,0.), length=5, width=5, height=7, centering="center bottom")
floor = Rectangle(position=(0.,0.,0.), length=10, width=50, height=30, centering="center top")
p_vel = player.velocity
np_norm = np.linalg.norm


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
	screen.draw_axes(10)


	#code
	player.update(dT, keys)
	dist = (screen.camera_position[2] - player.position[2])*6

	screen.draw_shape(floor, mode="lines", colour=(255-dist, 255-dist, 255-dist))
	screen.draw_shape(player, mode="lines", colour=(255-dist, 255-dist, 255-dist))
	
	# screen.draw_shape(floor, mode="fill")
	# print(player.velocity)

	screen.follow(player, offset=np.array([0,10,0]))
	screen.camera_position[2] = max(3,abs(np_norm(p_vel)))*3 + 20

	

	#tick end
	screen.update()

	if keys[pg.K_ESCAPE]:
		run = False
		print(player.position)
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break




