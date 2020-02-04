import pygame as pg
import pygame.locals as pglc
import numpy as np
from math import cos, sin, sqrt
import math, time
from scipy.spatial.distance import euclidean
from time import perf_counter
import modules.colour_maps as cmap
import modules.utils as utils



pg.init()


class PlotDrawer:
	#utility metods
	def __init__(self, size, plot_obj, axis_margin=(60, 100, 130, 100), title='', tick_length=15, **kwargs):
		self.plot_obj = plot_obj
		self._size = size
		self.axis_margin = axis_margin
		self.plotsurf = pg.surface.Surface(size)
		self.title = title
		self.tick_length = tick_length
		self.y_label = ''
		self.x_label = ''
		

	@property
	def size(self):
		return self._size


	@size.setter
	def size(self, val):
		self._size = val
		self.plotsurf = pg.surface.Surface(val)
		self._draw_plot()


	#drawing methods:
	def _draw_line(self, surf, p1, p2, width=1, colour=(255,255,255), aa=True):
		if aa: pg.draw.aaline(surf, colour, p1, p2)
		else: pg.draw.line(surf, colour, p1, p2, width)


	def _draw_lines(self, surf, p, width=1, colour=(255,255,255), aa=True):
		if aa: pg.draw.aalines(surf, colour, False, p)
		else: pg.draw.lines(surf, colour, False, p, width)

	def _draw_symbols(self, surf, ps, colour=(255,255,255), symbol='circle', size=4):
		if symbol == 'circle': 
			for p in ps:
				pg.draw.circle(surf, colour, (int(p[0]), int(p[1])), int(size))

	def _draw_symbol(self, surf, p, colour=(255,255,255), symbol='circle', size=4):
		if symbol == 'circle': 
			pg.draw.circle(surf, colour, (int(p[0]), int(p[1])), int(size))


	def _draw_text(self, surf, text, pos=(0,0), align='center', font='arial', size=20, rotation=0):
		text = str(text)

		f = pg.font.Font(pg.font.match_font(font), size)
		text_surf = f.render(text, True, (255,255,255))
		text_surf = pg.transform.rotate(text_surf, rotation)

		pos = np.asarray(pos)
		if align == 'center':
			pos[0] -= text_surf.get_size()[0]/2
			pos[1] -= text_surf.get_size()[1]/2
		if align == 'right center':
			pos[0] -= 0
			pos[1] -= text_surf.get_size()[1]/2
		if align == 'left center':
			pos[0] -= text_surf.get_size()[0]
			pos[1] -= text_surf.get_size()[1]/2

		surf.blit(text_surf, pos)


	def _draw_plot(self):
		'''
		Method that draws a plot.
		'''

		#get the offset in the left upper corner
		s = self.size
		offset = o = np.array([s[0]-self.axis_margin[0], s[1]-self.axis_margin[1]])
		opposite_offset = p = np.array([s[0] - s[0]+self.axis_margin[2], s[1] - s[1]+self.axis_margin[3]])

		#draw the square surrounding the plot
		self._draw_line(self.plotsurf, o, (p[0], o[1]))
		self._draw_line(self.plotsurf, o, (o[0], p[1]))
		self._draw_line(self.plotsurf, (p[0], o[1]), p)
		self._draw_line(self.plotsurf, (o[0], p[1]), p)

		#plot the data
		x_ran, x_ran_marg = (self.plot_obj.axes[0][1], self.plot_obj.axes[0][0]), (self.plot_obj.axes[2][1], self.plot_obj.axes[2][0])
		y_ran, y_ran_marg = self.plot_obj.axes[1], self.plot_obj.axes[3]

			#transform range function
		trans_x = lambda x: ((x - x_ran_marg[0])*(p[0]-o[0])/(x_ran_marg[1]-x_ran_marg[0])) + o[0]
		trans_y = lambda y: ((y - y_ran_marg[0])*(p[1]-o[1])/(y_ran_marg[1]-y_ran_marg[0])) + o[1]

		for series_x, series_y, colour, style in zip(self.plot_obj.x, self.plot_obj.y, self.plot_obj.colours.values(), self.plot_obj.plot_styles):
			x = series_x.copy()
			y = series_y.copy()

			x = trans_x(x)
			y = trans_y(y)

			if style == 'line':
				self._draw_lines(self.plotsurf, list(zip(x, y)), colour=colour)

			elif style == 'scatter':
				self._draw_lines(self.plotsurf, list(zip(x, y)), colour=colour)
				self._draw_symbols(self.plotsurf, list(zip(x, y)), colour=colour)

		#draw title
		title_pos = (s[0]/2, self.axis_margin[1]/2)
		self._draw_text(self.plotsurf, self.title, pos=title_pos)

		#draw_axis ticks
		p1, p2 = (trans_x(x_ran[0]), o[1]), (trans_x(x_ran[0]), o[1] + self.tick_length)
		self._draw_line(self.plotsurf, p1, p2)
		self._draw_text(self.plotsurf, round(x_ran[0], 1), size=15, pos=(p1[0], p1[1] + 2*self.tick_length))

		p1, p2 = (trans_x(x_ran[1]), o[1]), (trans_x(x_ran[1]), o[1] + self.tick_length)
		self._draw_line(self.plotsurf, p1, p2)
		self._draw_text(self.plotsurf, round(x_ran[1], 1), size=15, pos=(p1[0], p1[1] + 2*self.tick_length))

		p1, p2 = (p[0], trans_y(y_ran[0])), (p[0] - self.tick_length, trans_y(y_ran[0]))
		self._draw_line(self.plotsurf, p1, p2)
		self._draw_text(self.plotsurf, round(y_ran[0], 1), size=15, pos=(p1[0] - 2*self.tick_length, p1[1]), align='left center')

		p1, p2 =  (p[0], trans_y(y_ran[1])), (p[0] - self.tick_length, trans_y(y_ran[1]))
		self._draw_line(self.plotsurf, p1, p2)
		self._draw_text(self.plotsurf, round(y_ran[1], 1), size=15, pos=(p1[0] - 2*self.tick_length, p1[1]), align='left center')

		#draw axis labels
		p = ((p[0]-o[0])/2+o[0], s[1] - self.axis_margin[1]/2)
		self._draw_text(self.plotsurf, self.x_label, p)

		p = (self.axis_margin[1]/2, s[1]/2)
		self._draw_text(self.plotsurf, self.y_label, p, rotation=90)

	#show methods:
	def blit(self, surf_from, surf_to=None, origin=(0,0)):
		if surf_to == None:
			surf_to = self.disp

		surf_to.blit(surf_from, origin)


	def show(self):
		keepon = True
		disp = pg.display.set_mode(self.size, pglc.RESIZABLE)
		self.blit(self.plotsurf, disp)

		while keepon:
			ev = pg.event.get()
			keys = pg.key.get_pressed()

			if keys[pg.K_ESCAPE]:
				keepon = False
			for e in ev:
				if e.type == pg.QUIT:
					keepon = False
				if e.type == pg.VIDEORESIZE:
					self.size = e.dict['size']
					disp = pg.display.set_mode(self.size, pglc.RESIZABLE)
					self.blit(self.plotsurf, disp)

			pg.display.update()
					  
		pg.quit()





class Plot:
	'''
	Class that handles and contains data series.
	'''

	def __init__(self, colour_map=cmap.Rainbow(), size=(400,400)):
		self.x = []
		self.y = []
		self.label = []
		self.colours = {}
		self.id = []
		self.colour_map = colour_map
		self.plot_styles = []
		self.drawer = PlotDrawer(size, self)


	def plot(self, x, y, label=None, colour=None, style='line'):
		'''
		Method used to add data to the object.

		x, y - arrays of data, must be same dimension
		label - string with name of the series
		colour - tuple representing R, G, B
		plot_style - 'line' or 'scatter'
		'''

		x = np.asarray(x).flatten()
		y = np.asarray(y).flatten()

		assert(x.size == y.size)

		self.id.append(len(self.id))
		self.x.append(x)
		self.y.append(y)
		self.label.append(label)
		self.plot_styles.append(style)


		if not colour == None:
			self.colours[self.id[-1]] = colour

		self.axes = self._get_axes()


	def show(self):
		'''
		Draws and displays the plot using Renderer

		size - tuple specifying resolution of the renderer
		'''

		#set the colours:
		l = len(self.id)
		c = len(self.colours)
		for i in range(l):
			if not i in self.colours.keys():
				self.colours[i] = self.colour_map[i/(l-c)]

		
		self.drawer._draw_plot()
		self.drawer.show()


	def save(self, file):
		'''
		Method that saves the current plot.
		'''

		...


	def _get_axes(self, margin=0.05):
		'''
		Method that determines the ranges of the axes and the tick marks.
		'''

		#x-axis
		x_max = max([max(x) for x in self.x])
		x_min = min([min(x) for x in self.x])
		x_del = x_max - x_min

		x_max_marg = x_max + x_del * margin
		x_min_marg = x_min - x_del * margin


		#y-axis
		y_max = max([max(y) for y in self.y])
		y_min = min([min(y) for y in self.y])
		y_del = y_max - y_min

		y_max_marg = y_max + y_del * margin
		y_min_marg = y_min - y_del * margin

		return (x_min, x_max), (y_min, y_max), (x_min_marg, x_max_marg), (y_min_marg, y_max_marg)


