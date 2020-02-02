import numpy as np
from math import exp, log, sqrt, cos, sin, pi

class ForceField:
	###parameters:
	valence_bond = {
		'H_': 0.354,
		'H_b': 0.460,
		'He4+4': 0.849,
		'Li': 1.336,
		'Be3+2': 1.074,
		'B_3': 0.838,
		'B_2': 0.828,
		'C_3': 0.757,
		'C_R': 0.729,
		'C_2': 0.732,
		'C_1': 0.706,
		'N_3': 0.700,
		'N_R': 0.699,
		'N_2': 0.685,
		'N_1': 0.656,
		'O_3': 0.658,
		'O_3_z': 0.528,
		'O_R': 0.680,
		'O_2': 0.634,
		'O_1': 0.639,
		'F_': 0.668,
		'Ne4+4': 0.920,
	}

	valence_angle = {
		'H_': 180,
		'H_b': 83.5,
		'He4+4': 90,
		'Li': 180,
		'Be3+2': 109.47,
		'B_3': 109.47,
		'B_2': 120,
		'C_3': 109.47,
		'C_R': 120,
		'C_2': 120,
		'C_1': 180,
		'N_3': 106.7,
		'N_R': 120,
		'N_2': 111.2,
		'N_1': 180,
		'O_3': 104.51,
		'O_3_z': 146,
		'O_R': 110,
		'O_2': 120,
		'O_1': 180,
		'F_': 180,
		'Ne4+4': 90,
	}

	nonbond_distance = {
		'H_': 2.886,
		'H_b': 2.886,
		'He4+4': 2.362,
		'Li': 2.451,
		'Be3+2': 2.745,
		'B_3': 4.083,
		'B_2': 4.083,
		'C_3': 3.851,
		'C_R': 3.851,
		'C_2': 3.851,
		'C_1': 3.851,
		'N_3': 3.660,
		'N_R': 3.660,
		'N_2': 3.660,
		'N_1': 3.660,
		'O_3': 3.500,
		'O_3_z': 3.500,
		'O_R': 3.500,
		'O_2': 3.500,
		'O_1': 3.500,
		'F_': 3.364,
		'Ne4+4': 3.243,
	}

	nonbond_energy = {
		'H_': 0.044,
		'H_b': 0.044,
		'He4+4': 0.056,
		'Li': 0.025,
		'Be3+2': 0.085,
		'B_3': 0.180,
		'B_2': 0.180,
		'C_3': 0.105,
		'C_R': 0.105,
		'C_2': 0.105,
		'C_1': 0.105,
		'N_3': 0.069,
		'N_R': 0.069,
		'N_2': 0.069,
		'N_1': 0.069,
		'O_3': 0.060,
		'O_3_z': 0.060,
		'O_R': 0.060,
		'O_2': 0.060,
		'O_1': 0.060,
		'F_': 0.050,
		'Ne4+4': 0.042,
	}

	nonbond_scale = {
		'H_': 12,
		'H_b': 12,
		'He4+4': 15.24,
		'Li': 12,
		'Be3+2': 12,
		'B_3': 12.052,
		'B_2': 12.052,
		'C_3': 12.73,
		'C_R': 12.73,
		'C_2': 12.73,
		'C_1': 12.73,
		'N_3': 13.407,
		'N_R': 13.407,
		'N_2': 13.407,
		'N_1': 13.407,
		'O_3': 14.085,
		'O_3_z': 14.085,
		'O_R': 14.085,
		'O_2': 14.085,
		'O_1': 14.085,
		'F_': 14.762,
		'Ne4+4': 15.440,
	}

	effective_charge = {
		'H_': 0.712,
		'H_b': 0.712,
		'He4+4': 0.098,
		'Li': 1.026,
		'Be3+2': 1.565,
		'B_3': 1.755,
		'B_2': 1.755,
		'C_3': 1.912,
		'C_R': 1.912,
		'C_2': 1.912,
		'C_1': 1.912,
		'N_3': 2.544,
		'N_R': 2.544,
		'N_2': 2.544,
		'N_1': 2.544,
		'O_3': 2.300,
		'O_3_z': 2.300,
		'O_R': 2.300,
		'O_2': 2.300,
		'O_1': 2.300,
		'F_': 1.735,
		'Ne4+4': 0.194,
	}

	torsional_barrier_params = {
		'C_3': 2.119,
		'N_3': 0.450,
		'O_3': 0.018,
		'Si3': 1.225,
		'P_3': 2.400,
		'S_3': 0.484,
		'Ge3': 0.701,
		'As3': 1.500,
		'Se3': 0.335,
		'Sn3': 0.199,
		'Sb3': 1.100,
		'Te3': 0.300,
		'Pb3': 0.100,
		'Bi3': 1.000,
		'Po3': 0.300,
	}


	def get_atom_type(self, atom):
		el = atom.symbol 
		el += '_'*(len(el)==1)
		if atom.hybridisation != 0:
			el += str(atom.hybridisation)
		return el

	def get_energy(self, molecule):
		atoms = molecule.atoms

		#### E = E_r + E_theta + E_phi + E_ohm + E_vdm + E_el
		## E_r:
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
			print(e1, e2)
			#dist between a1 and a2
			r = a1.distance_to(a2)
			#bo between a1 and a2
			n = a1.bond_orders[a2]

			#components of eq dist
			rBO = -0.1332 * (r1+r2) * log(n)
			rEN = r1 * r2 * (sqrt(chi1)-sqrt(chi2))**2/(chi1*r1+chi2*r2)
			r12 = r1 + r2 + rBO + rEN

			
			#force constant
			k12 = 664.12 * Z1 * Z2 / r12**3 

			print(.5 * k12 * (r - r12)**2 * 4.2)

			#energy
			E_r += .5 * k12 * (r - r12)**2

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
			K123 = beta * Z1*Z3/r13**5 * r12*r23 * (r12*r23*(1-cos(t0)**2) - r13**2*cos(t0))

			C2 = 1/(4*sin(t0)**2)
			C1 = -4*C2*cos(t0)
			C0 = C2*(2*cos(t0)**2+1)

			E_theta += K123 * (C0 + C1*cos(theta) + C2*cos(2*theta))

		#E_phi:
		E_phi = 0
		for a1, a2, a3, a4, phi in molecule.get_unique_torsion_angles():
			e1 = self.get_atom_type(a1)
			e2 = self.get_atom_type(a2)
			e3 = self.get_atom_type(a3)
			e4 = self.get_atom_type(a4)

			Vbarr = 1
			if a2.hybridisation == a3.hybridisation == 3:
				V2, V3 = self.torsional_barrier_params[e2], self.torsional_barrier_params[e3]
				Vbarr = sqrt(V2*V3)
				n = 3
				phi0 = 180, 60 #two different phi0 possible

			if (a2.hybridisation == 3 and a3.hybridisation == 2) or (a2.hybridisation == 2 and a3.hybridisation == 3):
				Vbarr = 1
				n = 6
				phi0 = 0, 0

			if a2.hybridisation == a3.hybridisation == 2:
				U2 = (..., 2, 1.25, 0.7, 0.2, 0.1)[a2.period] # period starts at 1
				U3 = (..., 2, 1.25, 0.7, 0.2, 0.1)[a3.period]
				Vbarr = 5*sqrt(U2*U3)*(1+4.18*log(a2.bond_orders[a3]))
				n = 2
				phi0 = 180, 60

			#since two differen phi0 are possible, calculate both and return highest energy
			E_phi1 = 0.5*Vbarr * (1-cos(n*phi0[0])*cos(n*phi))
			E_phi2 = 0.5*Vbarr * (1-cos(n*phi0[1])*cos(n*phi))

			E_phi += max(E_phi1, E_phi2)


		print(f'TOTAL BOND STRETCHING ENERGY = {round(E_r*4.2,3)} kJ/mol')
		print(f'TOTAL ANGLE BENDING ENERGY = {round(E_theta*4.2,3)} kJ/mol')
		print(f'TOTAL TORSIONAL ENERGY = {round(E_phi*4.2,3)} kJ/mol')
		return E_r + E_theta + E_phi