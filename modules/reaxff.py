import numpy as np
from math import exp, log

class ForceField:
	def __init__(self, molecule):
		self.molecule = molecule

		self.general_params = [
		50.0,15.61,5.02,18.32,8.32,-8.90,1.94,
		-3.47,5.79,12.38,1.49,1.28,6.30,2.72,
		33.87,6.70,1.06,2.04,36.0,7.98,0.40,
		4.00,3.17,10.00,0.90,-1.14,2.17,1.69]

		self.atom_params = {
		'bond radii': {
			'C': [1.399, 1.266, 1.236],
			'H': [0.656]
			},
		'under/over coordination': {
			'C': [52.2, 29.4],
			'H': [117.5]
			},
		'Coulomb parameters': {
			'C': [7.41, 4.12, 0.69],
			'H': [9.14, 2.26, 0.37]
			},
		'heat increments': {
			'C': [218.6],
			'H': [54.3]
			}
		}

		self.bond_params = {
		'CC': [145.2,  0.318,  0.65, -0.097, 6.38, -0.26, 9.37, -0.391, 16.87],
		'CH': [183.8, -0.454, 12.80, -0.013, 7.65],
		'HH': [168.4, -0.310, 10.25, -0.016, 5.98]
		}


		self.set_bond_orders()

	def set_bond_orders(self):
		'''
		Method that creates a bond order matrix where element (i,j) represents
		the bond order between atoms i and j in self.molecule.atoms.
		'''

		atoms = self.molecule.atoms 
		bo = np.empty((len(atoms), len(atoms))) #create empty matrix

		#First set the uncorrected bond orders
		for i, a1 in enumerate(atoms):
			for j, a2 in enumerate(atoms):
				dist = a1.distance_to(a2) #get distance between atoms
				p = self.bond_params[''.join(sorted((a1.symbol, a2.symbol)))] #get the correct bond_params

				#case where a1 and a2 are both carbon
				if not a1.symbol == a2.symbol == 'C':
					#we only need the params for the sigma bonds (index 0):
					ra = (self.atom_params['bond radii'][a1.symbol][0] + self.atom_params['bond radii'][a2.symbol][0])/2
					bo[i,j] = exp(p[3] * (dist/ra)**p[4])

				if a1.symbol == a2.symbol == 'C':
					ra = (self.atom_params['bond radii'][a1.symbol][0] + self.atom_params['bond radii'][a2.symbol][0])/2
					rb = (self.atom_params['bond radii'][a1.symbol][1] + self.atom_params['bond radii'][a2.symbol][1])/2
					rc = (self.atom_params['bond radii'][a1.symbol][2] + self.atom_params['bond radii'][a2.symbol][2])/2
					bo[i,j] = exp(p[3]*(dist/ra)**p[4]) + exp(p[5]*(dist/rb)**p[6]) + exp(p[7]*(dist/rc)**p[8])

				if a1 == a2:
					bo[i,j] = 0

		#calculate the deviations between sum of valence and bond orders:
		bo_sums = np.sum(bo, axis=0)
		d = bo_sums - np.asarray([(1,4)[a.symbol == 'C'] for a in atoms])
		print(bo_sums)
		#correcting bond orders...
		l = self.general_params
		for i, a1 in enumerate(atoms):
			vali = (1,4)[a1.symbol == 'C']

			for j, a2 in enumerate(atoms):
				valj = (1,4)[a2.symbol == 'C']

				f5 = 1/(1 + exp(-l[2] * (l[3] * bo[i,j]**2 - d[i]) + l[4]))
				f4 = 1/(1 + exp(-l[2] * (l[3] * bo[i,j]**2 - d[j]) + l[4]))
				f3 = 1/l[1] * log(.5 * (exp(-l[1]*d[i]) + exp(-l[1]*d[j])))
				f2 = exp(-l[0]*d[i]) + exp(-l[0]*d[j])
				f1 = .5 * ((vali + f2)/(vali + f2 + f3) + (valj + f2)/(valj + f2 + f3))

				bo[i,j] = bo[i,j] * f1 * f4 * f5
		self.bond_orders = bo

