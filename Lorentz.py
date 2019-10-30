from modules.screen import *
from modules.shape import *
import math
import numpy as np
import random


#game setup
WIDTH, HEIGHT = SIZE = (600, 600)
disp = Screen3D(SIZE, camera_position=(0, 0, 600), camera_display_surface=(WIDTH/2, HEIGHT/2, 600))
clock = pg.time.Clock()
FPS = 100
run = True

#lorentz parameters
curr_pos = (100,-100,100)
sigma = 10
ro = 28
beta = 8/3

#main loop
updt = 0
while run:
	#tick prep
	updt += 1
	dT = clock.tick_busy_loop(FPS)/1000
	# disp.clear()

	#code
	disp.draw_axes(10)
	prev_pos = curr_pos


	x = dT * (sigma * (prev_pos[1] - prev_pos[0]))
	y = dT * (prev_pos[0] * (ro - prev_pos[2]) - prev_pos[1])
	z = dT * (prev_pos[0] * prev_pos[1] - beta * prev_pos[2])
	curr_pos = (x, y, z)
	print(dT)

	disp.draw_line((curr_pos, prev_pos))

	#tick end
	disp.update()
	for event in pg.event.get():
		if event.type == pg.QUIT:
			final = False