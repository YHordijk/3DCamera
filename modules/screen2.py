import pygame as pg
import numpy as np
from math import cos, sin
import math

class Screen3D:
	def __init__(self, size, camera_position=np.array([0.,0.,30.]), camera_orientation=(0.,0.,0.), project_type="perspective", bkgr_colour=(0,0,0)):
		self.camera_position = np.asarray(camera_position)
		self.camera_orientation = np.asarray(camera_orientation)
		self.project_type = project_type
		self.bkgr_colour = bkgr_colour
		self.size = self.width, self.height = size
		self.disp = pg.display.set_mode(self.size)
		self.disp.fill(self.bkgr_colour)

	def project(self, coord):
		a = coord.T
		c = self.camera_position.T
		t = self.camera_orientation
		e = np.array([self.width/2, self.height/2, 600])

		x_rot_mat = np.array([[1,0,0], [0, cos(t[0]), sin(t[0])], [0, -sin(t[0]), cos(t[0])]])
		y_rot_mat = np.array([[cos(t[1]), 0, -sin(t[1])], [0,1,0], [sin(t[1]), 0, cos(t[1])]])
		z_rot_mat = np.array([[cos(t[2]), sin(t[2]), 0], [-sin(t[2]), cos(t[2]), 0], [0,0,1]])
		z_rot_mat = np.array([[1,0,0], [0,1,0], [0,0,1]])

		d = x_rot_mat @ y_rot_mat @ z_rot_mat @ (a - c)

		f = np.array([[1, 0, e[0]/e[2]], [0, 1, e[1]/e[2]], [0, 0, 1/e[2]]]) @ d

		return int(round(f[0]/f[2])), int(round(f[1]/f[2]))

	def draw_pixel(self, pos, colour=(255,255,255)):
		try:
			pos = self.project(pos)
			self.disp.set_at(pos, colour)
		except:
			pass

	def draw_lines(self, poss, colour=(255,255,255), closed=True):
		try:
			poss = [self.project(pos) for pos in poss]
			pg.draw.aalines(self.disp, colour, closed, poss)
		except:
			raise

	def draw_line(self, poss, colour=(255,255,255)):
		try:
			poss = [self.project(pos) for pos in poss]
			pg.draw.aaline(self.disp, colour, poss[0], poss[1])
		except:
			raise

	def update(self):
		pg.display.update()

	def draw_shape(self, shape, colour=(255,255,255)):
		self.draw_lines(shape.points, colour)

	def draw_axes(self, length):
		self.draw_line([np.asarray((0,0,0)),np.asarray((length,0,0))], colour=(255,0,0))
		self.draw_line([np.asarray((0,0,0)),np.asarray((0,length,0))], colour=(0,255,0))
		self.draw_line([np.asarray((0,0,0)),np.asarray((0,0,length))], colour=(0,0,255))

	def clear(self):
		self.disp.fill(self.bkgr_colour)

	def look_at(self, point=None, obj=None):
		if point is None:
			point = obj.position

		x, y, z = point

		rotx = math.atan2(y, z)
		if z >= 0:
			roty = -math.atan2(x * math.cos(rotx), z)
		else:
			roty = math.atan2(x * math.cos(rotx), -z)
		rotz = math.atan2( math.cos(rotx), math.sin(rotx) * math.sin(roty))
		self.camera_orientation = np.array([roty, rotx, rotz])
		return self.camera_orientation
