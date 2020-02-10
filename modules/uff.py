import numpy as np
from math import exp, log, sqrt, cos, sin, pi
import os
import modules.utils as utils

class ForceField:
	def __init__(self):
		self.load_params()


	def load_params(self):
		file = os.getcwd() + r'\modules\data\UFF.prm'

		self.valence_bond = {}
		self.valence_angle = {}
		self.nonbond_distance = {}
		self.nonbond_scale = {}
		self.nonbond_energy = {}
		self.effective_charge = {}
		self.sp3_torsional_barrier_params = {}
		self.sp2_torsional_barrier_params = {}
		self.electro_negativity = {}

		with open(file, 'r') as f:
			for line in f.readlines():
				parts = line.split()
				if len(parts) > 0:
					if 'param' == parts[0]:
						e = parts[1]
						self.valence_bond[e] = float(parts[2])
						self.valence_angle[e] = float(parts[3])
						self.nonbond_distance[e] = float(parts[4])
						self.nonbond_energy[e] = float(parts[5])
						self.nonbond_scale[e] = float(parts[6])
						self.effective_charge[e] = float(parts[7])
						self.sp3_torsional_barrier_params[e] = float(parts[8])
						self.sp2_torsional_barrier_params[e] = float(parts[9])
						self.electro_negativity[e] = float(parts[10])


	def get_atom_type(self, atom):
		el = atom.symbol 
		el += '_'*(len(el)==1)
		if atom.ring == 'AR':
			el += 'R'
		elif atom.hybridisation != 0:
			el += str(atom.hybridisation)

		return el


	def get_energy(self, molecule, morse_potential=True, verbosity=0):
		atoms = molecule.atoms
		
		if verbosity > 1:
			utils.message('ATOM TYPES')
			utils.message(f'IDX | TYPE  | RING')
			for i, a in enumerate(atoms):
				utils.message(f'{i: <3} | {self.get_atom_type(a): <5} | {a.ring}')
				


		#### E = E_r + E_theta + E_phi + E_ohm + E_vdm + E_el
		## E_r:
		if verbosity > 2:
			utils.message('BOND STRETCH ENERGY')
			utils.message('ATOM 1 | ATOM 2 | BO  | BOND LEN | IDEAL LEN | ENERGY')

		E_r = 0
		for a1, a2 in molecule.get_unique_bonds():
			#some parameters for a1
			e1 = self.get_atom_type(a1)
			r1 = self.valence_bond[e1]
			chi1 = a1.electro_negativity
			Z1 = self.effective_charge[e1]

			#some parameters for a2
			e2 = self.get_atom_type(a2)
			r2 = self.valence_bond[e2]
			chi2 = a2.electro_negativity
			Z2 = self.effective_charge[e2]

			#dist between a1 and a2
			r = a1.distance_to(a2)
			#bo between a1 and a2
			n = a1.bond_orders[a2]
			if a1.ring == a2.ring == 'AR':
				n = 1.5

			if morse_potential:
				# r12 = r1 + r2 - (r1*r2*(sqrt(chi1) - sqrt(chi2))**2)/(chi1*r1 + chi2*r2)
				r12 = r1 + r2 - 0.1332 * (r1 + r2) * log(n) + r1*r2*(sqrt(chi1)-sqrt(chi2))**2/(chi1*r1+chi2*r2)
				k12 = 664.12 * Z1 * Z2 / r12**3 
				D12 = 70*n
				alpha = sqrt(k12/(2*D12))
				E_ri = D12*(exp(-alpha*(r-r12)) - 1)**2
			else:
				rBO = -0.1332 * (r1+r2) * log(n)
				rEN = r1 * r2 * (sqrt(chi1)-sqrt(chi2))**2/(chi1*r1+chi2*r2)
				r12 = r1 + r2 + rBO + rEN
				k12 = 664.12 * Z1 * Z2 / r12**3 
				E_ri = .5 * k12 * (r - r12)**2

			E_r += E_ri

			if verbosity > 2:
				utils.message(f'{e1: <6} | {e2: <6} | {n: <3.1f} | {r: <8.3f} | {r12: <8.3f} | {E_ri: <.3f}')


		## E_theta:
		

		E_theta = 0
		for a1, a2, a3, theta in molecule.get_unique_bond_angles():
			e1 = self.get_atom_type(a1)
			e2 = self.get_atom_type(a2)
			e3 = self.get_atom_type(a3)
			Z1 = self.effective_charge[e1]
			Z2 = self.effective_charge[e2]
			Z3 = self.effective_charge[e3]

			r12, r23, r13 = a1.distance_to(a2), a2.distance_to(a3), a1.distance_to(a3)

			beta = 664.12/(r12*r23)
			t0 = self.valence_angle[e2] * pi / 180		
			K123 = beta * Z1*Z3/r13**5 * r12*r23 * (3*r12*r23*(1-cos(t0)**2) - r13**2*cos(t0))
			C2 = 1/(4*sin(t0)**2)
			C1 = -4*C2*cos(t0)
			C0 = C2*(2*cos(t0)**2+1)
			E_theta_i = K123 * (C0 + C1*cos(theta) + C2*cos(2*theta))
			E_theta += E_theta_i

			

		#E_phi:
		if verbosity > 2:
			utils.message('BOND STRETCH ENERGY')
			utils.message('ATOM 1 | ATOM 2 | ATOM 3 | ATOM 4 | TORSION | IDEAL TORS | ENERGY')

		E_phi = 0
		for a1, a2, a3, a4, phi in molecule.get_unique_torsion_angles():
			e1 = self.get_atom_type(a1)
			e2 = self.get_atom_type(a2)
			e3 = self.get_atom_type(a3)
			e4 = self.get_atom_type(a4)

			Vbarr = 1
			if a2.hybridisation == a3.hybridisation == 3:
				V2, V3 = self.sp3_torsional_barrier_params[e2], self.sp3_torsional_barrier_params[e3]
				Vbarr = sqrt(V2*V3)
				n = 3
				phi0 = pi, pi/3 #two different phi0 possible

			if (a2.hybridisation == 3 and a3.hybridisation == 2) or (a2.hybridisation == 2 and a3.hybridisation == 3):
				Vbarr = 1
				n = 6
				phi0 = 0, 0

			if a2.hybridisation == a3.hybridisation == 2:
				U2 = self.sp2_torsional_barrier_params[e2] # period starts at 1
				U3 = self.sp2_torsional_barrier_params[e3]
				Vbarr = 5*sqrt(U2*U3)*(1+4.18*log(a2.bond_orders[a3]))
				n = 2
				phi0 = pi, 1.047198

			#since two differen phi0 are possible, calculate both and return highest energy
			E_phi1 = 0.5*Vbarr * (1-cos(n*phi0[0])*cos(n*phi))
			E_phi2 = 0.5*Vbarr * (1-cos(n*phi0[1])*cos(n*phi))
			E_phi_i = min(E_phi1, E_phi2)
			E_phi += E_phi_i

			if verbosity > 2:
				torsion = divmod(phi, phi0[0])[1] if phi0[0] is not 0 else 0
				utils.message(f'{e1: <6} | {e2: <6} | {e3: <6} | {e4: <6} | {torsion: <7.3f} | {phi0[0]: <10.3f} | {E_phi_i: <.3f}')

		#E_vdw:
		if verbosity > 2:
			utils.message('VAN DER WAALS ENERGY')
			utils.message('ATOM 1 | ATOM 2 | BOND LEN | FORCE CONST | ENERGY')

		E_vdw = 0
		for a1, a2 in molecule.get_unique_atom_pairs(): #cutoff is at 2 angstrom
			if a1.bond_dist_to(a2) > 2:
				e1 = self.get_atom_type(a1)
				e2 = self.get_atom_type(a2)

				D1 = self.nonbond_energy[e1]
				D2 = self.nonbond_energy[e2]
				D12 = sqrt(D1*D2)
				
				x = a1.distance_to(a2)
				x1 = self.nonbond_distance[e1]
				x2 = self.nonbond_distance[e2]
				x12 = .5*(x1+x2)

				E_vdw_i = D12 * (-2*(x12/x)**6 + (x12/x)**12)
				E_vdw += E_vdw_i

				if verbosity > 2:
					utils.message(f'{e1: <6} | {e2: <6} | {x: <8.3f} | {D12: <11.3f} | {E_vdw_i: <.3f}')




		if verbosity>0:
			utils.message(f'TOTAL BOND STRETCHING ENERGY = {round(E_r*4.2,3)} kJ/mol')
			utils.message(f'TOTAL ANGLE BENDING ENERGY = {round(E_theta*4.2,3)} kJ/mol')
			utils.message(f'TOTAL TORSIONAL ENERGY = {round(E_phi*4.2,3)} kJ/mol')
			utils.message(f'TOTAL VAN DER WAAL\'S ENERGY = {round(E_vdw*4.2,3)} kJ/mol')
			utils.message(f'TOTAL ENERGY = {round((E_r + E_theta + E_phi + E_vdw)*4.2,3)} kJ/mol')

		return E_r + E_theta + E_phi + E_vdw