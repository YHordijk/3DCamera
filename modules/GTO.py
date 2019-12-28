import math

class GTAO:
	def __init__(self, ao_type='1s', n=5):
		self.ao_type = ao_type
		self.coords = 0

		self.n = n
		self.l = list(range(n+1))[::-2]
		self.k = [k//2 for k in list(range(n+1))[::2]]

	def get_func(self):
		if self.ao_type == '1s':
			return lambda x, y, z: 0.08724638*math.exp(-0.151623*(x+y+z)) + 0.27181242724*math.exp(-0.851819*(x+y+z))
		return 


gtao = GTAO()
orb_1s = gtao.get_func()
print(orb_1s(0,0,0))


	# def exponent_term(self):
		# return np.exp(-)