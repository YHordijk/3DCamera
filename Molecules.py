import argparse
import modules.screen4 as scr 
import modules.molecule6 as mol 
import modules.basisset6 as bs
import math, os
import numpy as np
from time import perf_counter
import pubchempy as pcp
import pygame as pg

# parser = argparse.ArgumentParser(description='Molecule visualizer and molecular orbital calculator.')
# parser.add_argument('--molecule', help='Path to molecule or name of molecule. If name ends with ".pcp" it will be searched for on pubchem, else it will look for ".xyz" files.')
# parser.add_argument('-bs', '--basis_set', metavar='', type=str, default='STO-6G', help='Flag specifying the basis-set to be used.')
# parser.add_argument('-p', '--pre_render_dens', metavar='', type=bool, default=False, help='Flag specifying whether the molecular orbitals should be pre-rendered.')
# parser.add_argument('-R', '--resolution', metavar='', type=tuple, default=(1200,720), help='Flag specifying the width and height of screen.')
# parser.add_argument('-bc', '--background_colour', metavar='', type=tuple, default=(0, 0, 0), help='Flag specifying the red, green and blue components of the background-colour,')

# args = parser.parse_args()


pg.init()



#screen setup
# WIDTH, HEIGHT = SIZE = args.resolution
WIDTH, HEIGHT = SIZE = (1200,720)
# screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=args.background_colour)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=(0,0,0))
mols = [mol.Molecule('phenol.pcp', basis_set_type='STO-6G')]
# mols = [mol.Molecule(args.molecule, basis_set_type=args.basis_set)]
mol = mols[0]
# mol.add_hydrogens()
atoms = mol.atoms

bs.extended_huckel(mol) 
mos = mol.molecular_orbitals

# if args.pre_render_dens: screen.pre_render_densities(mos, points=10000)
screen.pre_render_densities(mos, points=10000)

pg.display.set_caption(', '.join([m.name.capitalize() for m in mols]))
clock = pg.time.Clock()
FPS = 120

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0

selected_atoms = set()

rot = np.array([0.,0.,0.])
zoom = 0

pg.key.set_repeat()
run = True

mo_numb = 0

while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT

	keys = pg.key.get_pressed()
	ev = pg.event.get()

	screen.clear()



	screen.draw_density(mos[mo_numb%len(mos)], 10000)

	[screen.draw_shape(m, wireframe=False, draw_atoms=True, draw_bonds=True, draw_hydrogens=True) for m in mols]
	[m.rotate(rot) for m in mols]

	screen.draw_axes(1)

	

	if keys[pg.K_ESCAPE]:
		run = False

	for e in ev:
		if e.type == pg.QUIT:
			run = False

		if e.type == pg.KEYDOWN:
			if e.key == pg.K_RIGHT:
				mo_numb += 1
			if e.key == pg.K_LEFT:
				mo_numb -= 1

		elif e.type == pg.MOUSEBUTTONDOWN:
			if e.button == 1:
				if not(keys[pg.K_LCTRL] or keys[pg.K_RCTRL]):
					pos = pg.mouse.get_pos()
					selected_atom = screen.get_atom_at_pos(pos)
					if not selected_atom == None:
						if selected_atom.selected:
							selected_atoms.discard(selected_atom)
							selected_atom.draw_colour = selected_atom.colour 
							selected_atom.selected = False
						else:
							selected_atoms.add(selected_atom)
							selected_atom.selected = True
							selected_atom.draw_colour = (247,240,111)
					else:
						selected_atoms = set()
						for mol in mols:
							mol.reset_colours()
							for a in mol.atoms:
								a.selected = False
			elif e.button == 4:
				if screen.camera_position[2] > 3:
					zoom = -dT * screen.camera_position[2] * 3
				else:
					screen.camera_position[2] = 3
			elif e.button == 5:
				if screen.camera_position[2] < 40:
					zoom = dT * screen.camera_position[2] * 3
				else:
					screen.camera_position[2] = 40
					

	rot *= 0.8
	move = pg.mouse.get_rel()
	if keys[pg.K_LCTRL] or keys[pg.K_RCTRL]:
		if pg.mouse.get_pressed()[2]:
				screen.camera_position[0] += move[0]/50
				screen.camera_position[1] += move[1]/50

		if pg.mouse.get_pressed()[0]:
			rot = np.asarray([move[1]/150, -move[0]/150, 0])


	screen.camera_position[2] += zoom
	zoom *= 0.8


	#tick end
	l = len(selected_atoms)
	if l > 0:
		sa = list(selected_atoms)
		if l == 1:
			atom = list(selected_atoms)[0]

			if len(mols) > 1:
				for mol in mols:
					if atom in mol.atoms:
						index = mol.atoms.index(atom)
			else:
				index = mols[0].atoms.index(atom)

			screen.display_text(f'  {atom.symbol}{index}  ', (10,10))

		elif l == 2:
			atoms = list(selected_atoms)
			index = []
			for atom in atoms:
				if len(mols) > 1:
					for mol in mols:
						if atom in mol.atoms:
							index.append(atoms.index(atom))
				else:
					index.append(mols[0].atoms.index(atom))

			screen.display_text(f'  {atoms[0].symbol}{index[0]}, {atoms[1].symbol}{index[1]}: {round(atoms[0].distance_to(atoms[1]), 3)} (A)  ', (10,10))

		elif l == 3:
			atoms = list(selected_atoms)
			index = []
			for atom in atoms:
				if len(mols) > 1:
					for mol in mols:
						if atom in mol.atoms:
							index.append(atoms.index(atom))
				else:
					index.append(mols[0].atoms.index(atom))

			a1, a2, a3 = atoms
			if a2 in a1.bonds and a3 in a1.bonds:
				screen.display_text(f'  {atoms[1].symbol}{index[1]}, {atoms[0].symbol}{index[0]}, {atoms[2].symbol}{index[2]}: {round(mols[0].bond_angle(a2, a1, a3, in_degrees=True), 1)} (deg)  ', (10,10))
			elif a1 in a2.bonds and a3 in a2.bonds:
				screen.display_text(f'  {atoms[0].symbol}{index[0]}, {atoms[1].symbol}{index[1]}, {atoms[2].symbol}{index[2]}: {round(mols[0].bond_angle(a1, a2, a3, in_degrees=True), 1)} (deg)  ', (10,10))
			elif a1 in a3.bonds and a2 in a3.bonds:
				screen.display_text(f'  {atoms[0].symbol}{index[0]}, {atoms[2].symbol}{index[2]}, {atoms[1].symbol}{index[1]}: {round(mols[0].bond_angle(a1, a3, a2, in_degrees=True), 1)} (deg)  ', (10,10))
			else:
				screen.display_text(f'  {atoms[0].symbol}{index[0]}, {atoms[1].symbol}{index[1]}, {atoms[2].symbol}{index[2]}  ', (10,10))


	screen.update()
