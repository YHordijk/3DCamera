import pygame as pg
import numpy as np
from math import cos, sin, pi
import math


class Screen2D:
	def __init__(self, SIZE, rangex=(-1,1), rangey=(-1,1), COLOR=(0,0,0)):
		self.WIDTH, self.HEIGHT = self.SIZE = SIZE
		self.COLOR = COLOR
		self.rangex = rangex
		#rangey must be inverted because the y-axis in pygame is upside down
		self.rangey = -rangey[1], -rangey[0]
		self.disp = pg.display.set_mode(self.SIZE)
		self.disp.fill(self.COLOR)

	def transform(self, pos):
		x, y = pos[0], -pos[1]
		x = (x - self.rangex[0])*self.WIDTH/(self.rangex[1]-self.rangex[0])
		y = (y - self.rangey[0])*self.HEIGHT/(self.rangey[1]-self.rangey[0])
		return int(x), int(y)

	def check_pixel(self, pos):
		try:
			# print(pos)
			return self.disp.get_at(pos)
		except IndexError as e:
			return False

	def draw_pixel(self, pos, COLOR=(255,255,255), ALPHA=None):
		pos = self.transform(pos) 
		if not ALPHA is None:
			#target  colors
			ta, tb, tc = COLOR
			#current colors of pixel
			print(self.check_pixel(pos))
			ca, cb, cc = self.check_pixel(pos)[:3]
			#background color
			ba, bb, bc = self.COLOR
			#difference between target and background
			da, db, dc = ta-ba, tb-bb, tc-bc
			#get new color
			na, nb, nc = (ca + da*(ALPHA/100)), (cb + db*(ALPHA/100)), (cc + dc*(ALPHA/100))
			#check if legitimate color assignment
			na = min(255, na, ba)
			na = max(0, na, ba)
			nb = min(255, nb, bb)
			nb = max(0, nb, bb)
			nc = min(255, nc, bc)
			nc = max(0, nc, bc)
			COLOR = (na, nb, nc)
	
		#transform coordinates to new basis (rangex, rangey)
		self.disp.set_at(pos, COLOR)

	def update(self, final=False):
		pg.display.flip()

		while final:
			for event in pg.event.get():
				if event.type == pg.QUIT:
					final = False

	def clear(self):
		self.disp.fill(self.COLOR)

class Screen3D:
	def __init__(self, SIZE, camera_position=(0,0,0), camera_orientation=(0,0,0), camera_display_surface=None, camera_ortho_properties=(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1), project_type="perspective", COLOR=(0,0,0)):
		self.camera_position = camera_position
		self.camera_orientation = camera_orientation
		self.project_type = project_type
		self.camera_display_surface = camera_display_surface if not camera_display_surface is None else (SIZE[0], SIZE[1], 100)
		self.camera_ortho_properties = camera_ortho_properties
		self.WIDTH, self.HEIGHT = self.SIZE = SIZE
		self.COLOR = COLOR
		self.disp = pg.display.set_mode(self.SIZE)
		self.disp.fill(self.COLOR)

	def project(self, point):
		if self.project_type == "perspective":
			ax, ay, az = point
			cx, cy, cz = self.camera_position
			tx, ty, tz = self.camera_orientation
			ex, ey, ez = self.camera_display_surface
			
			x, y, z = ax-cx, ay-cy, az-cz
			sx, sy, sz = sin(tx), sin(ty), sin(tz)
			cx, cy, cz = cos(tx), cos(ty), cos(tz)

			dx = cy * (sz*y + cz*x) - sy*z
			dy = sx * (cy*z + sy*(sz*y + cz*x)) + cx * (cz*y - sz*x)
			dz = cx * (cy*z + sy*(sz*y + cz*x)) - sx * (cz*y - sz*x)

			bx = ez*dx/dz + ex
			by = ez*dy/dz + ey

			return (int(bx), int(by))

		if self.project_type == "orthographic":
			p = self.camera_ortho_properties
			ax, ay, az = point
			point = np.array(([ax, ay, az, 1]))
			P = np.array(([ 2/(p[1]-p[0]), 0, 0, (-p[1]-p[0])/(p[1]-p[0])],
						  [0,  2/(p[3]-p[2]), 0, (-p[3]-p[2])/(p[3]-p[2])],
						  [0, 0, -2/(p[5]-p[4]), (-p[5]-p[4])/(p[5]-p[4])],
						  [0, 0, 0, 1]))
			B = P @ np.asarray(point)
			return (int(B[0]), int(B[1]))


	def draw_pixel(self, point, color=(255,255,255)):
		pos = self.project(point)
		self.disp.set_at(pos, color)

	def draw_line(self, points, color=(255,255,255)):
		pos = [self.project(p) for p in points]
		pg.draw.lines(self.disp, color, True, pos)

	def draw_shape(self, shape, color=(255,255,255)):
		pos = [self.project(p) for p in shape.points]
		pg.draw.lines(self.disp, color, True, pos)

	def update(self):
		pg.display.flip()

	def clear(self):
		self.disp.fill(self.COLOR)

	def draw_axes(self, length):
		self.draw_line([(0,0,0),(length,0,0)], color=(255,0,0))
		self.draw_line([(0,0,0),(0,length,0)], color=(0,255,0))
		self.draw_line([(0,0,0),(0,0,length)], color=(0,0,255))

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
