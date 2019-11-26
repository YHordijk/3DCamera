import pygame as pg
import numpy as np
from math import cos, sin, pi


class Shape:
	def __init__(self, rotation=(0.,0.,0.), position=(0.,0.,0.)):
		self.rotation = np.asarray(rotation)
		self.position = np.asarray(position)

	def convert_points_to_array(self, points):
		if type(points) is list:
			new_points = []
			for p in self.points:
				new_points.append(np.array((p[0], p[1], p[2])))
			points = new_points
		return points

	def convert_points_to_list(self, points):
		if type(points) is np.ndarray:
			points.tolist()


	def rotate(self, points, rotation):
		# self.points = self.convert_points_to_array(self.points)

		r = rotation[0]
		Rx = np.array(([	  1, 	  0,	   0],
					   [	  0, cos(r), -sin(r)],
					   [      0, sin(r),  cos(r)]))

		r = rotation[1]
		Ry = np.array(([ cos(r),  	   0, sin(r)],
					   [ 	  0, 	   1,	   0],
					   [-sin(r), 	   0, cos(r)]))

		r = rotation[2]
		Rz = np.array(([ cos(r), -sin(r), 	   0],
					   [ sin(r),  cos(r), 	   0],
					   [ 	  0, 	   0, 	   1]))

		new_points = [Rx @ Ry @ Rz @ p for p in points]

		return new_points



class Rectangle(Shape):
	def __init__(self, height=10, width=10, length=10, centering="center", **kwargs):
		super().__init__(**kwargs)
		self.height = height
		self.width = width
		self.length = length
		self.centering = centering
		self.closed = False
		self.type = 'flat'

	@property
	def points(self):
		h = self.height
		w = self.width
		l = self.length

		if self.centering == "center":
			offset = np.array([w/2, h/2, l/2])
		elif self.centering == "center top":
			offset = np.array([w/2, h, l/2])
		elif self.centering == "center bottom":
			offset = np.array([w/2, 0, l/2])
		else:
			offset = 0

		points = np.array(([0,0,0], 
						   [0,h,0], 
						   [0,h,l], 
						   [0,0,l], 
						   [0,0,0], 
						   [w,0,0], 
						   [w,h,0], 
						   [0,h,0], 
						   [w,h,0], 
						   [w,h,l], 
						   [w,0,l], 
						   [w,0,0], 
						   [w,0,l], 
						   [0,0,l], 
						   [0,h,l], 
						   [w,h,l])) - offset

		return self.rotate(points, self.rotation) + np.array(([self.position]))

	def faces(self, double_sided=True):
		h = self.height
		w = self.width
		l = self.length

		if self.centering == "center":
			offset = np.array([w/2, h/2, l/2])
		elif self.centering == "center top":
			offset = np.array([w/2, h, l/2])
		else:
			offset = 0

		faces = np.array(([[0,0,0],[w,0,0],[w,0,l],[0,0,l]],
						  [[w,0,l],[w,h,l],[0,h,l],[0,0,l]],
						  [[0,h,0],[w,h,0],[w,h,l],[0,h,l]])) - offset

		return [self.rotate(face, self.rotation) + np.array(([self.position])) for face in faces]



class Square(Shape):
	def __init__(self, height=10, width=10, centering="center", **kwargs):
		super().__init__(**kwargs)
		self.height = height
		self.width = width
		self.centering = centering
		self.closed = True
		self.type = 'flat'

	@property
	def points(self):
		h = self.height
		w = self.width

		if self.centering == "center":
			offset = np.array([w/2, h/2, 0])
		else:
			offset = 0

		points = np.array(([0,0,0], 
						   [0,h,0], 
						   [w,h,0], 
						   [w,0,0])) - offset

		return self.rotate(points, self.rotation) + np.array(([self.position]))

	def faces(self, double_sided=True):
		h = self.height
		w = self.width

		if self.centering == "center":
			offset = np.array([w/2, h/2, 0])
		elif self.centering == "center top":
			offset = np.array([w/2, h, 0])
		else:
			offset = 0

		faces = np.array(([[[0,0,0],[w,0,0],[w,h,0],[0,h,0]],])) - offset
		faces = [self.rotate(face, self.rotation) + np.array(([self.position])) for face in faces]
		return faces


class Cube(Shape):
	def __init__(self, edge, centering="center", **kwargs):
		super().__init__(**kwargs)
		self.edge = edge
		self.centering = centering
		self.closed = False
		self.type = 'flat'

	@property
	def points(self):
		e = self.edge

		if self.centering == "center":
			offset = e/2
		else:
			offset = 0
		
		points = np.array(([0,0,0], 
						   [0,e,0], 
						   [0,e,e], 
						   [0,0,e], 
						   [0,0,0], 
						   [e,0,0], 
						   [e,e,0], 
						   [0,e,0], 
						   [e,e,0], 
						   [e,e,e], 
						   [e,0,e], 
						   [e,0,0], 
						   [e,0,e], 
						   [0,0,e], 
						   [0,e,e], 
						   [e,e,e])) - offset

		points = self.rotate(points, self.rotation) + np.array(([self.position]))
		return points

	@property
	def verteces(self):
		e = self.edge

		if self.centering == "center":
			offset = e/2
		else:
			offset = 0

		points = np.array(([0,0,0], 
						   [0,e,0], 
						   [e,e,0], 
						   [e,0,0], 
						   [0,0,e], 
						   [0,e,e], 
						   [e,e,e], 
						   [e,0,e],)) - offset

		return self.rotate(points, self.rotation) + np.array(([self.position]))


class Atom(Shape):
		def __init__(self, element='H', **kwargs):
			super().__init__(**kwargs)
			self.element = element
			self.dictionaries = {
			'C':{'colour': (34, 34, 34), 'radius': 67},
			'H':{'colour': (255, 255, 255), 'radius': 53},
			'O':{'colour': (255, 22, 0), 'radius': 48},
			'N':{'colour': (22, 33, 255), 'radius': 56},
			'S':{'colour': (225, 225, 48), 'radius': 100},
			'Ca':{'colour': (61, 255, 0), 'radius': 180},
			}

			self.colour = self.dictionaries[self.element]['colour']
			self.radius = self.dictionaries[self.element]['radius']

class Molecule(Shape):
	def __init__(self, atoms=[], **kwargs):
		super().__init__(**kwargs)
		self.atoms = atoms
		self.type = 'molecule'
		self.scale_factor = 4

	def load_xyz(self, xyz):
		elements = np.loadtxt(xyz, skiprows=2, usecols=0, dtype=str)
		coords = np.loadtxt(xyz, skiprows=2, usecols=(1,2,3), dtype=float)

		self.atoms = []
		for e, c in zip(elements, coords):
			self.atoms.append(Atom(element=e, position=c))

		for a in self.atoms:
			a.radius = int(a.radius * self.scale_factor)

	def load_sdf(self, sdf):
		pass

	def update_atoms(self):
		points = [a.position for a in self.atoms]
		points = self.rotate(points, self.rotation)

		for a, p in zip(self.atoms, points):
			a.position = p