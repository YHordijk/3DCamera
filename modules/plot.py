import pygame as pg
import pygame.locals as pglc
import numpy as np
from math import cos, sin, sqrt
import math, time
from scipy.spatial.distance import euclidean
from time import perf_counter
import modules.colour_maps as cmap
import modules.utils as utils

## MARCHING SQUARES DATA
'''
	
   0____0___ 1
	|		|
	|		|
   3|		| 1
	|		|
	|_______|
   3	2	 2

s = a, b, c or d
vertex_lt[i] = (s11, s21, s12, s22)
'''

vertex_lt = [(-1, -1, -1, -1),
			 ( 3,  2, -1, -1),
			 ( 2,  1, -1, -1),
			 ( 3,  1, -1, -1),
			 ( 0,  1, -1, -1),
			 ( 0,  3,  1,  2),
			 ( 0,  2, -1, -1),
			 ( 0,  3, -1, -1),
			 ( 0,  3, -1, -1),
			 ( 0,  2, -1, -1),
			 ( 3,  2,  0,  1),
			 ( 0,  1, -1, -1),
			 ( 3,  1, -1, -1),
			 ( 2,  1, -1, -1),
			 ( 3,  2, -1, -1),
			 (-1, -1, -1, -1)]

side_offset_lt = np.array([(0,0), (1,0), (0,1), (0,0)])





class Series:
	'''
	Class holding data for data series
	'''

	def __init__(self, x, y, label='', style='line', z=None, colour=None, colour_map=None, width=None):
		self.x = x
		self.y = y
		self.label = label
		self.style = style
		# if style is contour or heatmap we need the z values
		if style in ['contour', 'heatmap']:
			self.z = z
			self.colour_map = colour_map
		else:
			self.colour = colour
			self.width = width

class Plot:
	'''
	Class that handles and contains data series.
	'''

	def __init__(self, size=(400,400), colour_map=cmap.Rainbow(), axis_margin=(60, 100, 130, 100), title='', tick_length=15, axes=None, plot_margin=0.05):
		self._size = size
		self.colour_map = colour_map
		self.axis_margin = axis_margin
		self.title = title
		self.tick_length = tick_length
		self.y_label = ''
		self.x_label = ''
		self.custom_axes = axes
		self.plot_margin = plot_margin

		self.plotsurf = pg.surface.Surface(size)
		self.series = []



	##### PLOTTING PORTION
	def plot_3d(self, x, y, z, label=None, colour_map=cmap.CoolWarm(), style='heatmap'):
		'''
		Method used to add data to the object.

		x, y - vectors with x and y data
		z - matrix holding the data in de third dimension

		label - string with name of the series
		colour - tuple representing R, G, B
		plot_style - 'line' or 'scatter'
		'''

		x = np.asarray(x).flatten()
		y = np.asarray(y).flatten()
		z = np.asarray(z)
		assert(z.size == x.size*y.size)
		assert(style in ['contour', 'heatmap'])

		self.series.append(Series(x, y, z=z, label=label, colour_map=colour_map, style=style))


	def plot(self, x, y, label=None, colour=None, style='line', width=1):
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
		assert(style in ['line', 'scatter'])

		self.series.append(Series(x, y, label=label, colour=colour, style=style, width=width))


	def save(self, file, resolution=None):
		'''
		Method that saves the current plot.
		'''

		if resolution == None:
			resolution = self.size
		...


	def _get_axes(self):
		'''
		Method that determines the ranges of the axes and the tick marks.
		'''

		if 'heatmap' in [s.style for s in self.series] or 'contour' in [s.style for s in self.series]: margin = 0
		else: margin = self.plot_margin

		#x-axis
		x_max = max([max(s.x) for s in self.series])
		x_min = min([min(s.x) for s in self.series])
		x_del = x_max - x_min

		x_max_marg = x_max + x_del * margin
		x_min_marg = x_min - x_del * margin


		#y-axis
		y_max = max([max(s.y) for s in self.series])
		y_min = min([min(s.y) for s in self.series])
		y_del = y_max - y_min

		y_max_marg = y_max + y_del * margin
		y_min_marg = y_min - y_del * margin

		return (x_min, x_max), (y_min, y_max), (x_min_marg, x_max_marg), (y_min_marg, y_max_marg)


	def clear(self):
		self.series = []

	##### DRAWING PORTION
	@property
	def size(self):
		return self._size


	@size.setter
	def size(self, val):
		self._size = val
		self.plotsurf = pg.surface.Surface(val)
		self._draw_plot()

	def _transform_to_plot_range(self, surf, pos):
		if self.custom_axes == None: 
			axes = self._get_axes()
		else: 
			axes = self.custom_axes

		s = surf.get_size()
		o = np.array([0,0])
		p = np.array(s)

		x_ran = axes[0]
		y_ran = axes[1]

		x, y = pos

		x = ((x - x_ran[0])*(s[0])/(x_ran[1]-x_ran[0]))
		y = ((y - y_ran[0])*(s[1])/(y_ran[1]-y_ran[0]))

		return x, y


	#drawing methods:
	def _draw_contours(self, surf, ser, thresh):
		'''
		Method that draws the contour lines of a heatmap or other 3d plot
		'''

		assert(ser.style in ['contour', 'heatmap'])

		lines = []

		
		x, y = ser.x.copy(), ser.y.copy()
		dx = x[1] - x[0]
		dy = y[1] - y[0]

		z = ser.z.copy()

		vo = side_offset_lt

		lx, ly = ser.x.size, ser.y.size

		for i in range(lx-1):
			for j in range(ly-1):
				#get the values at the corners of the cell
				vals = np.array([z[i,j], z[i+1,j], z[i+1,j+1], z[i,j+1]])
				#get which corners are above thresh
				verteces = np.where(vals>=thresh, 1, 0)
				#get the case (convert binary verteces to base 10)
				case = sum([n*2**i for i,n in enumerate(verteces[::-1])])
				#get the sides of cell to connect
				case_verteces = vertex_lt[case]

				c1 = case_verteces[:2]
				c2 = case_verteces[2:4]

				if not -1 in c1:
					c1a = c1[0]
					c1b = c1[1]

					# if c1a == 0:
					# 	if c1b == 1:
					# 		p1 = 



					p1 = vo[c1a] * np.array([dx,dy])
					p2 = vo[c1b] * np.array([dx,dy])

					# p1 *= np.array([dx,dy])
					# p2 *= np.array([dx,dy])

					p1 += np.array([x[i], y[j]])
					p2 += np.array([x[i], y[j]])

					p1 = self._transform_to_plot_range(surf, p1)
					p2 = self._transform_to_plot_range(surf, p2)
					
					self._draw_line(surf, p1, p2)

					if not -1 in c2:
						c1a = c2[0]
						c1b = c2[1]

						#start at the offsets
						p1 = vo[c1a] * np.array([dx,dy])
						p2 = vo[c1b] * np.array([dx,dy])

						p1 += np.array([x[i], y[i]])
						p2 += np.array([x[i], y[i]])

						p1 = self._transform_to_plot_range(surf, p1)
						p2 = self._transform_to_plot_range(surf, p2)

						
						self._draw_line(surf, p1, p2)


	def _draw_line(self, surf, p1, p2, width=1, colour=(255,255,255), aa=False):
		if aa: pg.draw.aaline(surf, colour, p1, p2)
		else: pg.draw.line(surf, colour, p1, p2, width)


	def _draw_lines(self, surf, p, width=1, colour=(255,255,255), aa=False):
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


	def _get_ui_surf(self):
		'''
		Method that draws a plot.
		'''

		ui_surf = pg.surface.Surface(self.size)


		if self.custom_axes == None: 
			axes = self._get_axes()
		else: 
			axes = self.custom_axes

		#get the offset in the left upper corner
		s = self.size
		offset = o = np.array([s[0]-self.axis_margin[0], s[1]-self.axis_margin[1]])
		opposite_offset = p = np.array([s[0] - s[0]+self.axis_margin[2], s[1] - s[1]+self.axis_margin[3]])

		#plot the data
		x_ran, x_ran_marg = (axes[0][1], axes[0][0]), (axes[2][1], axes[2][0])
		y_ran, y_ran_marg = axes[1], axes[3]

		#transform range function
		trans_x = lambda x: ((x - x_ran_marg[0])*(p[0]-o[0])/(x_ran_marg[1]-x_ran_marg[0])) + o[0]
		trans_y = lambda y: ((y - y_ran_marg[0])*(p[1]-o[1])/(y_ran_marg[1]-y_ran_marg[0])) + o[1]	
	
		#draw the square surrounding the plot
		self._draw_line(ui_surf, o, (p[0], o[1]))
		self._draw_line(ui_surf, o, (o[0], p[1]))
		self._draw_line(ui_surf, (p[0], o[1]), p)
		self._draw_line(ui_surf, (o[0], p[1]), p)			

		#draw title
		title_pos = (s[0]/2, self.axis_margin[1]/2)
		self._draw_text(ui_surf, self.title, pos=title_pos)

		#draw_axis ticks
		p1, p2 = (trans_x(x_ran[0]), o[1]), (trans_x(x_ran[0]), o[1] + self.tick_length)
		self._draw_line(ui_surf, p1, p2)
		self._draw_text(ui_surf, round(x_ran[0], 1), size=15, pos=(p1[0], p1[1] + 2*self.tick_length))

		p1, p2 = (trans_x(x_ran[1]), o[1]), (trans_x(x_ran[1]), o[1] + self.tick_length)
		self._draw_line(ui_surf, p1, p2)
		self._draw_text(ui_surf, round(x_ran[1], 1), size=15, pos=(p1[0], p1[1] + 2*self.tick_length))

		p1, p2 = (p[0], trans_y(y_ran[0])), (p[0] - self.tick_length, trans_y(y_ran[0]))
		self._draw_line(ui_surf, p1, p2)
		self._draw_text(ui_surf, round(y_ran[0], 1), size=15, pos=(p1[0] - 2*self.tick_length, p1[1]), align='left center')

		p1, p2 =  (p[0], trans_y(y_ran[1])), (p[0] - self.tick_length, trans_y(y_ran[1]))
		self._draw_line(ui_surf, p1, p2)
		self._draw_text(ui_surf, round(y_ran[1], 1), size=15, pos=(p1[0] - 2*self.tick_length, p1[1]), align='left center')

		#draw axis labels
		pos = ((p[0]-o[0])/2+o[0], s[1] - self.axis_margin[1]/2)
		self._draw_text(ui_surf, self.x_label, pos)

		pos = (self.axis_margin[1]/2, s[1]/2)
		self._draw_text(ui_surf, self.y_label, pos, rotation=90)

		#draw transparency rect with colour (0,255,0)
		pg.draw.rect(ui_surf, (0,255,0), pg.Rect(p+1, o-p-1))
		ui_surf.set_colorkey((0,255,0))

		return ui_surf


	def _get_plot_surf(self, size):
		p_surf = pg.surface.Surface(size)

		for ser in self.series:
			x = ser.x.copy()
			y = ser.y.copy()

			x, y = self._transform_to_plot_range(p_surf, (x,y))
			if ser.style in ['line', 'scatter']:
				self._draw_lines(p_surf, list(zip(x, y)), colour=ser.colour, width=ser.width)
				if ser.style == 'scatter':
					self._draw_symbols(p_surf, list(zip(x, y)), colour=ser.colour)

			elif ser.style in ['heatmap', 'contour']:
				z = ser.z.copy()

				hm_surf = pg.surface.Surface(z.shape)
				pg.surfarray.blit_array(hm_surf, ser.colour_map[z])

				hm_surf = pg.transform.scale(hm_surf, size)
				p_surf.blit(hm_surf, (0,0))

				if ser.style == 'contour':
					self._draw_contours(p_surf, ser, 0)

		return p_surf


	def _draw_plot(self):
		s = self.size
		offset = o = np.array([s[0]-self.axis_margin[0], s[1]-self.axis_margin[1]])
		opposite_offset = p = np.array([s[0] - s[0]+self.axis_margin[2], s[1] - s[1]+self.axis_margin[3]])

		plot_offset = ((o-p)*self.plot_margin).astype(int)
		plot_size = (o-p) - plot_offset*2
		plot_size = int(plot_size[0]), int(plot_size[1])

		p_surf = self._get_plot_surf(plot_size)
		ui_surf = self._get_ui_surf()


		p_surf = pg.transform.scale(p_surf, plot_size)
		p_surf = pg.transform.flip(p_surf, False, True)

		self.plotsurf.blit(p_surf, p+plot_offset)
		self.plotsurf.blit(ui_surf, (0,0))

		# self.plotsurf = p_surf



	#show methods:
	def blit(self, surf_from, surf_to=None, origin=(0,0)):
		if surf_to == None:
			surf_to = self.disp

		surf_to.blit(surf_from, origin)


	def show(self):
		pg.init()

		# data series part
		colourable_series = [s for s in self.series if s.style in ['line', 'scatter']]
		c = sum([s.colour is not None for s in colourable_series])
		l = len(self.series)
		for i, s in enumerate(colourable_series):
			if s.colour is None:
				s.colour = self.colour_map.get_colour(i/(l-c))

		self._draw_plot()

		#drawing part
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