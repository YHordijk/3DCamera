from math import exp, cos

class ForceField:
	def __init__(self):
		pass


	def get_energy(self, molecule):
		pass


class MM2(ForceField):
	def __init__(self):
		self.params = {
						'stretching': 
							{
							'CC': 
								{'ks':4.4,
								'l0':1.523},

							'CH': 
							 	{'ks':4.6,
							 	'l0':1.113}
							},
					    'bending':
							{
							'CCC': 
								[
									{'k0':0.45,
									't0':109.5,
									'k0cb':0.34},

									{'k0':0.45,
									't0':109.5,
									'k0cb':0.29},

									{'k0':0.45,
									't0':109.5,
									'k0cb':0.34}
								],

							'CCH': 
								[
									{'k0':0.36,
									't0':109.4},

									{'k0':0.36,
									't0':109.4},

									{'k0':0.36,
									't0':110.0}
								],

							'CHH': 
								[
									{'k0':0.32,
									't0':109.4},

									{'k0':0.32,
									't0':109.0}
								]
							},
						'stretch-bend': 
							{
							'CCC':
								{'ks':0.12},

							'CCH':
								{'ks':0.09},

							'CHH':
								{'ks':0.00}
							},
						'torsion':
							{
							'CCCC':
								{'V1':0.20,
								'V2':0.27,
								'V3':0.093},

							'CCCH':
								{'V1':0.00,
								'V2':0.00,
								'V3':0.267},

							'CCHH':
								{'V1':0.00,
								'V2':0.00,
								'V3':0.237}
							},
						'vanderwaals':
							{
							'CC':
								{'Sr*':3.80,
								'e':0.044},

							'CH':
								{'Sr*':3.34,
								'e':0.046},

							'HH':
								{'Sr*':3.00,
								'e':0.047}
							}
						}


	def get_energy(self, molecule):
		if not ('C' in molecule.get_elements() and 'H' in molecule.get_elements() and len(molecule.get_elements()) == 2):
			raise ValueError('Only C and H atoms are compatible with the MM2 forcefield.')

		params = self.params

		#Electrostatic interactions:
		Ees = 0
		for a1, a2 in molecule.get_unique_atom_pairs():
			Ees += (a1.nuc_charge + a2.nuc_charge)/a1.distance_to(a2)

		#stretch energy:
		Es = 0
		stretching_params = params['stretching']
		for a1, a2 in molecule.get_unique_bonds():
			p = stretching_params[''.join(sorted((a1.element,a2.element)))]
			l = a1.distance_to(a2)
			Es += 71.94 * p['ks'] * (l - p['l0'])**2 * (1 - 2.00 * (l - p['l0']))


		#bending energy:
		Eb = 0
		bending_params = params['bending']
		for a1, a2, a3, angle in molecule.get_unique_bond_angles(in_degrees=True):
			elements = a1.element+a2.element+a3.element
			elements = ''.join(sorted(elements))

			if elements == 'CCC':
				atype = len(molecule.get_bonded_atoms(a2, 'C')) - 2
			elif elements == 'CCH' or elements == 'CHH':
				atype = len(molecule.get_bonded_atoms(a2, 'C')) - 1

			# print(a1, a2, a3, atype+1, elements, len(molecule.get_bonded_atoms(a2, 'C')))

			p = bending_params[elements][atype]
			Eb += 0.021914 * p['k0'] * (angle - p['t0'])**2 * (1 + 7e-8 * (angle - p['t0'])**4)


		#stretch-bend energy
		Esb = 0
		stretch_bend_params = params['stretch-bend']
		for a1, a2, a3, angle in molecule.get_unique_bond_angles(in_degrees=True):
			elements = a1.element+a2.element+a3.element
			elements = ''.join(sorted(elements))

			if elements == 'CCC':
				atype = len(molecule.get_bonded_atoms(a2, 'C')) - 2
			elif elements == 'CCH' or elements == 'CHH':
				atype = len(molecule.get_bonded_atoms(a2, 'C')) - 1

			la = a2.distance_to(a1)
			lb = a2.distance_to(a3)

			psb = stretch_bend_params[elements]
			pb = bending_params[elements][atype]

			psa = stretching_params[''.join(sorted((a1.element,a2.element)))]
			psb = stretching_params[''.join(sorted((a3.element,a2.element)))]


			Esb += 2.51124 * psb['ks'] * (angle - pb['t0']) * ((la - psa['l0']) + (lb - psb['l0']))


		#vdw energy:
		Ev = 0
		vdw_params = params['vanderwaals']
		for a1, a2 in molecule.get_unique_atom_pairs():
			elements = ''.join(sorted((a1.element,a2.element)))
			p = vdw_params[elements]
			P =p['Sr*']/a1.distance_to(a2)
			if elements == 'CH':
				P /= 0.915
			Ev += p['e'] * (2.90e5 * exp(-12.50/(P)) - 2.25*P**6)


		#Torsion energy:
		Et = 0
		tors_params = params['torsion']
		for a1, a2, a3, a4, angle in molecule.get_unique_torsion_angles(in_degrees=True):
			elements = a1.element+a2.element+a3.element+a4.element
			elements = ''.join(sorted(elements))

			p = tors_params[elements]

			Et += p['V1']/2 * (1 + cos(angle)) + p['V2']/2 * (1 - cos(2*angle)) + p['V3']/2 * (1 + cos(3*angle))

		return Es + Ees + Eb + Esb + Ev + Et

