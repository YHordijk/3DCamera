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

class Cube(Shape):
	def __init__(self, edge, centering="center", **kwargs):
		super().__init__(**kwargs)
		self.edge = edge
		self.centering = centering

	@property
	def points(self):
		if self.centering == "center":
			offset = self.edge/2
		else:
			offset = 0
		e = self.edge
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

# c = Cube(edge=5, position=(0,0,0))
# print(c.points)
# c = Cube(edge=5, position=(0,0,10))
# print(c.points)

# print(type(c).__bases__)