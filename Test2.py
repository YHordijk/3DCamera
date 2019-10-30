from modules.screen import *
from modules.shape import *
import math
import numpy as np
import random


#game setup
WIDTH, HEIGHT = SIZE = (600, 600)
disp = Screen3D(SIZE, camera_position=(3, 30, 101), camera_orientation=(0,0,0), camera_display_surface=(WIDTH/2, HEIGHT/2, 600))
clock = pg.time.Clock()
FPS = 60
run = True

#main loop
updt = -100
while run:
	#tick prep
	updt += 1
	dT = clock.tick_busy_loop(FPS)/1000
	# disp.clear()

	#code
	disp.look_at((0,0,0))
	disp.draw_axes(10)
	
	# disp.camera_orientation=(0.2,updt/10/FPS,0)
	for i in range(-100,100):
		disp.draw_pixel((updt, (i+updt)**2/100, -i))

	#tick end
	disp.update()
	for event in pg.event.get():
		if event.type == pg.QUIT:
			final = False