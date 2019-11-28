import pygame as pg
import numpy as np
from math import cos, sin
import math
from scipy.spatial.distance import euclidean
from time import perf_counter


class Screen3D:
	def __init__(self, size, camera_position=np.array([0.,0.,30.]), camera_orientation=(0.,0.,0.), project_type="perspective", bkgr_colour=(0,0,0)):
		self.camera_position = np.asarray(camera_position)
		self.camera_orientation = np.asarray(camera_orientation)
		self.project_type = project_type
		self.bkgr_colour = bkgr_colour
		self.size = self.width, self.height = size
		self.disp = pg.display.set_mode(self.size)
		self.disp.fill(self.bkgr_colour)

		self.trail = []

	def follow(self, target, offset=0):
		delta = (target.position - self.camera_position + offset)/4
		self.camera_position += delta


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
		except Exception as e:
			pass

	def draw_pixels(self, poss, colour=(255,255,255)):
		draw = self.draw_pixel
		[draw(pos, colour) for pos in poss]

	def draw_lines(self, poss, colour=(255,255,255), closed=True, width=1):
		try:
			proj = self.project
			poss = [proj(pos) for pos in poss]
			pg.draw.aalines(self.disp, colour, closed, poss, width)

		except:
			raise

	def draw_single_bond(self, poss, colour=(255,255,255), width=1):
		try:
			h = self.height + 200
			w = self.width + 200
			poss = [self.project(pos) for pos in poss]
			if (-200 <= poss[0][1] <= h and -200 <= poss[0][0] <= w and -200 <= poss[1][1] <= h and -200 <= poss[1][0] <= w):
				# print(poss, euclidean(poss[0], poss[1]))
				pg.draw.line(self.disp, self.bkgr_colour, poss[0], poss[1], width+3)
				pg.draw.line(self.disp, colour, poss[0], poss[1], width)
		except:
			raise

	def draw_double_bond(self, poss, colour=(255,255,255), width=1):
		try:
			h = self.height + 200
			w = self.width + 200
			poss = np.asarray([self.project(pos) for pos in poss])
			if (-200 <= poss[0][1] <= h and -200 <= poss[0][0] <= w and -200 <= poss[1][1] <= h and -200 <= poss[1][0] <= w):
				d = width
				perp = poss - poss[0]
				
				perp = np.asarray([perp[0], (perp[1][1], -perp[1][0])])[1]
				perp = d * perp / np.linalg.norm(perp)

				pg.draw.line(self.disp, self.bkgr_colour, poss[0]-perp, poss[1]-perp, d+3)
				pg.draw.line(self.disp, self.bkgr_colour, poss[0]+perp, poss[1]+perp, width+3)

				pg.draw.line(self.disp, colour, poss[0]-perp, poss[1]-perp, d)
				pg.draw.line(self.disp, colour, poss[0]+perp, poss[1]+perp, d)
		except:
			pass

	def draw_triple_bond(self, poss, colour=(255,255,255), width=1):
		try:
			h = self.height + 200
			w = self.width + 200
			poss = np.asarray([self.project(pos) for pos in poss])
			if (-200 <= poss[0][1] <= h and -200 <= poss[0][0] <= w and -200 <= poss[1][1] <= h and -200 <= poss[1][0] <= w):
				d = width
				perp = poss - poss[0]
				
				perp = np.asarray([perp[0], (perp[1][1], -perp[1][0])])[1]
				perp = 1.5 * d * perp / np.linalg.norm(perp)

				pg.draw.line(self.disp, self.bkgr_colour, poss[0]-perp, poss[1]-perp, d+3)
				pg.draw.line(self.disp, self.bkgr_colour, poss[0], poss[1], d+3)
				pg.draw.line(self.disp, self.bkgr_colour, poss[0]+perp, poss[1]+perp, width+3)

				pg.draw.line(self.disp, colour, poss[0]-perp, poss[1]-perp, d)
				pg.draw.line(self.disp, colour, poss[0], poss[1], d)
				pg.draw.line(self.disp, colour, poss[0]+perp, poss[1]+perp, d)
		except:
			pass


	def draw_circle(self, center, radius, colour=(255,255,255), width=0):
		pos = self.project(center)
		pg.draw.circle(self.disp, colour, pos, radius, width)


	def draw_polygon(self, points, colour=(255,255,255)):
		try:
			proj = self.project
			dist = int((self.camera_position[2] - points[0][1])*2)
			points = [proj(point) for point in points]
			colour = (255-dist,255-dist,255-dist)
			pg.draw.polygon(self.disp, colour, points)

		except:
			raise

	def update(self):
		pg.display.update()


	def draw_shape(self, shape, colour=(255,255,255), double_sided=False, mode="fill", draw_atoms=True, draw_bonds=True, colour_bonds=True):
		if shape.type == 'flat':
			if mode == "fill":
				faces = shape.faces(double_sided=double_sided)
				[self.draw_polygon(face, colour) for face in faces]
			if mode == "lines":
				self.draw_lines(shape.points, colour, shape.closed)

		elif shape.type == 'molecule':
			cam_pos = self.camera_position
			original_coords = shape.coords
			dists = np.asarray([shape.distance(cam_pos, x) for x in shape.coords])
			deltas = np.asarray(shape.coords - self.camera_position)
			indices = np.argsort(dists)[::-1]
			coords = shape.coords[indices]
			original_atoms = shape.atom_types
			atoms = shape.atom_types[indices]
			dists = dists[indices]
			deltas = deltas[indices]

			radii = shape._atom_radii
			orders = shape.bond_orders
			shape_pos = shape.position
			colours = shape._atom_colours

			prev_indices = []
			for i, a, c, d, delt in zip(indices, atoms, coords, dists, deltas):
				prev_indices.append(i)
				if delt[2] < 0:
					bonds = shape.bonds[i].copy()
					order = orders[i].copy()

					connected_atoms = original_coords[bonds]

					if draw_bonds:
						if colour_bonds:
							d2 = lambda x, y: (x - y)/2
							for o, a2, c1 in zip(order, bonds, connected_atoms):
								if not a2 in prev_indices:
									if o != shape.bond_orders[i][bonds.index(a2)]: print(f'{i}: {o}, {a2}: {shape.bond_orders[i][bonds.index(a2)]}'); print(orders[i])
									if o == 1:
										self.draw_single_bond((c1+shape_pos-d2(c1,c), c1+shape_pos), width=int(75/d), colour=colours[original_atoms[a2]])	
										self.draw_single_bond((c+shape_pos, c+shape_pos+d2(c1,c)), width=int(75/d), colour=colours[a])

									if o == 2:
										self.draw_double_bond((c1+shape_pos-d2(c1,c), c1+shape_pos), width=int(75/d), colour=colours[original_atoms[a2]])
										self.draw_double_bond((c+shape_pos, c+shape_pos+d2(c1,c)), width=int(75/d), colour=colours[a])

									if o == 3:
										self.draw_triple_bond((c1+shape_pos-d2(c1,c), c1+shape_pos), width=int(75/d), colour=colours[original_atoms[a2]])
										self.draw_triple_bond((c+shape_pos, c+shape_pos+d2(c1,c)), width=int(75/d), colour=colours[a])
						else:	
							[self.draw_line((c+shape_pos, c1+shape_pos), width=int(75/d)+3, colour=self.bkgr_colour) for c1 in connected_atoms]
							[self.draw_line((c+shape_pos, c1+shape_pos), width=int(75/d)) for c1 in connected_atoms]

					if draw_atoms:
						self.draw_circle(c+shape_pos, int(radii[a]/d * shape.scale)+1, self.bkgr_colour, width=2)
						self.draw_circle(c+shape_pos, int(radii[a]/d * shape.scale), colours[a])			

	def draw_axes(self, length):
		self.draw_line([np.asarray((0,0,0)),np.asarray((length,0,0))], colour=(255,0,0))
		self.draw_line([np.asarray((0,0,0)),np.asarray((0,length,0))], colour=(0,255,0))
		self.draw_line([np.asarray((0,0,0)),np.asarray((0,0,length))], colour=(0,0,255))

	def clear(self):
		self.disp.fill(self.bkgr_colour)

	def move(self, vector):
		a = np.asarray(vector).T
		t = self.camera_orientation
		c = self.camera_position.T
		x_rot_mat = np.array([[1,0,0], [0, cos(t[0]), sin(t[0])], [0, -sin(t[0]), cos(t[0])]])
		y_rot_mat = np.array([[cos(t[1]), 0, -sin(t[1])], [0,1,0], [sin(t[1]), 0, cos(t[1])]])
		z_rot_mat = np.array([[cos(t[2]), sin(t[2]), 0], [-sin(t[2]), cos(t[2]), 0], [0,0,1]])
		z_rot_mat = np.array([[1,0,0], [0,1,0], [0,0,1]])
		result = x_rot_mat @ y_rot_mat @ z_rot_mat @ (a)
		self.camera_position += result.T

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
