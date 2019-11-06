import modules.shape as shape
import pygame as pg
import numpy as np

class Player(shape.Rectangle):
	def __init__(self, velocity=(0.,0.,0.), acceleration=(0., 0., 0.), acceleration_on_move=20, terminal_velocity=20, break_force=2, **kwargs):
		super().__init__(**kwargs)
		self.acceleration = np.asarray(acceleration)
		self.velocity = np.asarray(velocity)
		self.acceleration_on_move = acceleration_on_move
		self.terminal_velocity = terminal_velocity
		self.break_force = break_force

	def collision_detection(self, objs):
		pass

	def update(self, dT, keys=None):
		if keys is None:
			keys = pg.key.get_pressed()

		p = self.position
		v = self.velocity
		a = self.acceleration
		aom = self.acceleration_on_move

		moved = False
		if keys[pg.K_a]:
			a[0] = aom
			moved = True
		if keys[pg.K_d]:
			a[0] = -aom
			moved = True
		if keys[pg.K_w]:
			a[1] = aom
			moved = True
		if keys[pg.K_s]:
			a[1] = -aom
			moved = True
		if not moved:
			a = np.array([0.,0.,0.])
			v /= 1 + (self.break_force * dT)

		v += a * dT
		v = np.minimum(v, self.terminal_velocity)
		v = np.maximum(v, -self.terminal_velocity)
		p += v * dT

		self.position = p
		self.velocity = v
		self.acceleration = a


