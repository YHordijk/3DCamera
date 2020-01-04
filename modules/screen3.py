import pygame as pg
import numpy as np
from math import cos, sin, sqrt
import math
from scipy.spatial.distance import euclidean
from time import perf_counter
import modules.colour_maps as cmap


class Screen3D:
	def __init__(self, size, camera_position=np.array([0.,0.,30.]), camera_orientation=(0.,0.,0.), project_type="perspective", bkgr_colour=(0,0,0)):
		self.camera_position = np.asarray(camera_position)
		self.camera_orientation = np.asarray(camera_orientation)
		self.project_type = project_type
		self.bkgr_colour = bkgr_colour
		self.size = self.width, self.height = size
		self.disp = pg.display.set_mode(self.size)
		self.disp.fill(self.bkgr_colour)


	def follow(self, target, offset=0):
		delta = (target.position - self.camera_position + offset)/4
		self.camera_position += delta


	def display_text(self, text, pos):
		f = pg.font.Font(pg.font.get_default_font(), 20)
		surf = f.render(text, True, (255,255,255))
		self.disp.blit(surf, pos)


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

	def project_array(self, array):
		#accepts array in the form:
		'''
		array = [x0, y0, z0]
				[x1, y1, z1]
					...
				[xn, yn, zn]
		'''

		a = array
		c = self.camera_position.T
		t = self.camera_orientation
		e = np.array([self.width/2, self.height/2, 600])

		x_rot_mat = np.array([[1,0,0], [0, cos(t[0]), sin(t[0])], [0, -sin(t[0]), cos(t[0])]])
		y_rot_mat = np.array([[cos(t[1]), 0, -sin(t[1])], [0,1,0], [sin(t[1]), 0, cos(t[1])]])
		z_rot_mat = np.array([[cos(t[2]), sin(t[2]), 0], [-sin(t[2]), cos(t[2]), 0], [0,0,1]])
		z_rot_mat = np.array([[1,0,0], [0,1,0], [0,0,1]])


		d = x_rot_mat @ y_rot_mat @ z_rot_mat @ (a - c).T

		f = np.array([[1, 0, e[0]/e[2]], [0, 1, e[1]/e[2]], [0, 0, 1/e[2]]]) @ d
		return np.vstack((f[0]/f[2], f[1]/f[2])).T.astype(int)

	def draw_pixel(self, pos, colour=(255,255,255)):
		try:
			pos = self.project(pos)
			self.disp.set_at(pos, colour)
		except Exception as e:
			pass
	
	def draw_pixels(self, poss, colour=(255,255,255), colour_func=None, colour_array=None):
		set_at = self.disp.set_at
		if type(colour_array) is np.ndarray:
			if type(poss) is np.ndarray:
				poss = self.project_array(poss)
				[set_at(pos, clr) for pos, clr in zip(poss, colour_array)]
			else:
				draw = self.draw_pixel
				[draw(pos, colour_func(pos)) for pos in poss]

		elif colour_func != None:
			if type(poss) is np.ndarray:
				poss = self.project_array(poss)
				[self.disp.set_at(pos, colour_func(pos)) for pos in poss]
			else:
				draw = self.draw_pixel
				[draw(pos, colour_func(pos)) for pos in poss]

		else:
			if type(poss) is np.ndarray:
				poss = self.project_array(poss)
				[set_at(pos, colour) for pos in poss]
			else:
				draw = self.draw_pixel
				[draw(pos, colour_func(pos)) for pos in poss]



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


	def get_atom_at_pos(self, pos):
		pos = np.asarray(pos)
		dists = sorted([(euclidean(pos, p), a) for a, p in self.atom_draw_pos.items()], key=lambda x: x[0])
		atom = None
		for dist, a in dists:
			if dist < self.atom_draw_rad[a]:
				try:
					if a.distance_to(self.camera_position) < atom.distance_to(self.camera_position):
						atom = a
				except:
					atom = a
		return atom


	def draw_density(self, molecule, points=50000, colour_map=cmap.BlackWhite()):
		if not hasattr(self, '_dens_pos'):
			print(f'Calculating density of {molecule.name} with {molecule.basis_set.basis_type} ...')
			samples = 50*points
			rang = np.amax([np.abs(atom.coords) for atom in molecule.atoms]) + 4

			x, y, z = ((np.random.randint(-rang*10000, rang*10000, size=samples)/10000), (np.random.randint(-rang*10000, rang*10000, size=samples)/10000), (np.random.randint(-rang*10000, rang*10000, size=samples)/10000))
			d = molecule.get_orb_density(np.asarray((x, y, z)).T).flatten()
			index = d.argsort()[::-1]
			colours = colour_map[d].T

			# index = np.arange(0, samples)
			# index = np.where(abs(d) > np.amax(d)/2, index, 0)

			# index = index[index > 0]

			x, y, z, colours = x[index][0:points], y[index][0:points], z[index][0:points], colours[index][0:points]
			self._dens_pos, self._dens_colours = np.asarray((x, y, z)).T, colours

		self.draw_pixels(self._dens_pos, colour_array=self._dens_colours)

	def draw_shape(self, shape, colour=(255,255,255), double_sided=False, mode="fill", draw_atoms=True, draw_bonds=True, colour_bonds=True, draw_hydrogens=True, wireframe=False):
		if shape.type == 'flat':
			if mode == "fill":
				faces = shape.faces(double_sided=double_sided)
				[self.draw_polygon(face, colour) for face in faces]
			if mode == "lines":
				self.draw_lines(shape.points, colour, shape.closed)

		elif shape.type == 'molecule':
			atom_draw_pos = {}
			atom_draw_rad = {}
			p = shape.position

			cam_pos = self.camera_position #reference to camera position
			if draw_hydrogens:
				atoms = shape.atoms #reference to the atoms
			else:
				atoms = shape.get_by_element('H', blacklist=True)
			coords = np.asarray([a.coords for a in atoms]) #coordinates converted to np array

			dists = np.asarray([a.distance_to((cam_pos*2)) for a in atoms])#calculate dists to determine order of drawing
			indices = np.argsort(dists)[::-1] #determine order of drawing by sorting the dists and reversing
			atoms = np.asarray(atoms)[indices] #sort atoms by distance to cam_pos
			dists = dists[indices]
			deltas = coords - cam_pos #determine delta distance to determine the 
			deltas = deltas[indices]

			d2 = lambda x, y: (x - y)/2

			prev_indices = []
			for i, a1 in enumerate(atoms):
				if not wireframe:
					width = int(150/dists[i])

					prev_indices.append(a1)
					# if deltas[i][2] < 0:
					if True:
						c1 = a1.coords
						if draw_bonds:
							if colour_bonds:
								for a2 in a1.bonds:
									if not a2.element == 'H' or draw_hydrogens:
										if not a2 in prev_indices:
											c2 = a2.coords
											if a1.bond_orders[a2] == 1:
												self.draw_single_bond((c1 + p, c1 + p + d2(c2,c1)), width=width, colour=a1.draw_colour)
											elif a1.bond_orders[a2] == 2:
												self.draw_double_bond((c1 + p, c1 + p + d2(c2,c1)), width=width, colour=a1.draw_colour)
											elif a1.bond_orders[a2] == 3:
												self.draw_triple_bond((c1 + p, c1 + p + d2(c2,c1)), width=width, colour=a1.draw_colour)


						if draw_atoms:
							rad = int(a1.radius/dists[i] * shape.scale)
							self.draw_circle(c1+p, rad+1, self.bkgr_colour)
							self.draw_circle(c1+p, rad, a1.draw_colour)
							atom_draw_pos[a1] = self.project(c1+p)
							atom_draw_rad[a1] = rad*1.1

						if draw_bonds:
							if colour_bonds:
								for a2 in a1.bonds:
									if not a2.element == 'H' or draw_hydrogens:
										if not a2 in prev_indices:
											c2 = a2.coords
											if a1.bond_orders[a2] == 1:
												self.draw_single_bond((c2 + p - d2(c2,c1), c2 + p), width=width, colour=a2.draw_colour)	
											elif a1.bond_orders[a2] == 2:
												self.draw_double_bond((c2 + p - d2(c2,c1), c2 + p), width=width, colour=a2.draw_colour)	
											elif a1.bond_orders[a2] == 3:
												self.draw_triple_bond((c2 + p - d2(c2,c1), c2 + p), width=width, colour=a2.draw_colour)	
				elif wireframe:
					c1 = a1.coords
					for a2 in a1.bonds:
						if not a2 in prev_indices:
							c2 = a2.coords
							if a1.bond_orders[a2] == 1:
								self.draw_single_bond((c1 + p, c1 + p + d2(c2,c1)), width=3, colour=a1.draw_colour)
								self.draw_single_bond((c2 + p - d2(c2,c1), c2 + p), width=3, colour=a2.draw_colour)
							elif a1.bond_orders[a2] == 2:
								self.draw_double_bond((c1 + p, c1 + p + d2(c2,c1)), width=3, colour=a1.draw_colour)
								self.draw_double_bond((c2 + p - d2(c2,c1), c2 + p), width=3, colour=a2.draw_colour)	
							elif a1.bond_orders[a2] == 3:
								self.draw_triple_bond((c1 + p, c1 + p + d2(c2,c1)), width=3, colour=a1.draw_colour)
								self.draw_triple_bond((c2 + p - d2(c2,c1), c2 + p), width=3, colour=a2.draw_colour)


			self.atom_draw_pos = atom_draw_pos
			self.atom_draw_rad = atom_draw_rad

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
					print(shape)
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


	# def on_click(self, x, y):



	def draw_axes(self, length=1):
		self.draw_single_bond([np.asarray((0,0,0)),np.asarray((length,0,0))], colour=(255,0,0))
		self.draw_single_bond([np.asarray((0,0,0)),np.asarray((0,length,0))], colour=(0,255,0))
		self.draw_single_bond([np.asarray((0,0,0)),np.asarray((0,0,length))], colour=(0,0,255))

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
