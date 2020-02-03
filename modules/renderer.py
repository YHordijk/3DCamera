import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"
import pygame as pg
import numpy as np
from math import cos, sin, pi, tan
import math, sys
try:
	import modules.colour_maps as cmap
except:
	import colour_maps as cmap
import time
import decimal as dec


class Renderer:
	def __init__(self, resolution=(0,0), rangex=None, rangey=None, colour_map=cmap.CoolWarm()):
		self.resolution = resolution
		self.rangex = rangex
		self.rangey = rangey
		self.colour_map = colour_map
		self.pixel_array = np.zeros(resolution)


	def clear(self):
		self.pixel_array = np.zeros(resolution)


	def map_rgb_to_surface(self, array):
		return cmap.Ocean[array]


	@property
	def resolution(self):
		return self._resolution


	@resolution.setter
	def resolution(self, val):
		self._resolution = val
		self.disp = pg.surface.Surface(val)
	

	def blit_array(self, array):
		pg.surfarray.blit_array(self.disp, self.colour_map[self.pixel_array])


	def show(self, clickable=False, range_set=False):
		self.blit_array(self.pixel_array)
		dest = pg.display.set_mode(self.resolution)
		dest.blit(self.disp, (0,0))
		pg.display.flip()
		keepon = True
		prev_press = 0
		prev_pos = []
		drag_surf = pg.surface.Surface(self.resolution)
		drag_surf.set_colorkey((255,0,0))
		resize = False

		while keepon:
			if clickable:
				if pg.mouse.get_pressed()[0] and not prev_press:
					x, y = pg.mouse.get_pos()
				
					if range_set:
						if len(prev_pos) == 1:
							prev_pos.append((x,y))
							tx = lambda x: x*(self.rangex[1]-self.rangex[0]) / self.resolution[0] + self.rangex[0]
							ty = lambda y: y*(self.rangey[1]-self.rangey[0]) / self.resolution[1] + self.rangey[0]

							self.rangex = dec.Decimal(tx(min(prev_pos[0][0], prev_pos[1][0]))), dec.Decimal(tx(max(prev_pos[0][0], prev_pos[1][0])))
							self.rangey = ty(min(prev_pos[0][1], prev_pos[1][1])), ty(max(prev_pos[0][1], prev_pos[1][1]))

							keepon = False
							resize = True

						if len(prev_pos) == 0:
							prev_pos.append((x,y))

					else:
						x = x*(self.rangex[1]-self.rangex[0]) / self.resolution[0] + self.rangex[0]
						y = y*(self.rangey[1]-self.rangey[0]) / self.resolution[1] + self.rangey[0]
						print(f' x = {(x):.3e},  y = {(y):.3e}')

					time.sleep(0.01)

				prev_press = pg.mouse.get_pressed()[0]

				if len(prev_pos) == 1:
					x, y = pg.mouse.get_pos()
					drag_surf.fill((255,0,0))
					pos = [prev_pos[0], (prev_pos[0][0], y), (x,y), (x, prev_pos[0][1])]
					pg.draw.lines(drag_surf, (0,0,0), True, pos, 2)
					dest.blit(self.disp, (0,0))
					dest.blit(drag_surf, (0,0))

				pg.display.update()

			if pg.key.get_pressed()[pg.K_ESCAPE]:
				keepon = False
			for event in pg.event.get():
				if event.type == pg.QUIT:
					  keepon = False
		pg.quit()

		return resize


	def save(self, path):
		self.blit_array(self.pixel_array)
		pg.image.save(self.disp, path)


	def input_array(self, array):
		self.pixel_array = array
		self.resolution = array.shape


	def input_pos(self, poss, auto_size=True):
		x, y = np.hsplit(poss, 2)
		if auto_size:
			self.rangex = x.min(), x.max()
			self.rangey = y.min(), x.max()
		x, y = self.transform_to_disp((x,y))

		poss = np.append(x, y)

		pa = self.pixel_array
		for p in poss:
			pa[p[0],p[1]] += 1

		self.pixel_array = pa

		return pa


	def transform_to_disp(self, pos):
		if type(pos) is tuple:
			x, y = pos
		else:
			x, y = np.hsplit(pos,2)[0], np.hsplit(pos,2)[1]

		x = ((x - self.rangex[0])*(self.resolution[0]-1)/(self.rangex[1]-self.rangex[0]))
		y = ((y - self.rangey[0])*(self.resolution[1]-1)/(self.rangey[1]-self.rangey[0]))

		return x.astype(int), y.astype(int)





def draw_cmap_sample(colour_map):
	res = (600,200)
	array = np.empty(res)
	s = Renderer(res, colour_map=colour_map)
	for y in range(res[1]):
		for x in range(res[0]):
			array[x,y] = x
	s.input_array(array)
	s.show(clickable=False)

