import modules.screen3 as scr 
import modules.molecule3 as mol 
import modules.forcefields as ff
import math, os
import numpy as np
from time import perf_counter
import pubchempy as pcp
import pygame as pg

pg.init()

# mols = [mol.Molecule(molecule_file='ethane.pcp', warning_level=1, scale=400)]

mols = [mol.Molecule(molecule_file=os.getcwd() + f'\\Molecules\\ethane.xyz', warning_level=1, scale=400)]

# mols = [mol.Molecule(molecule_file='Glucose.pcp', warning_level=1, position=[5,0,0], scale=400),
# 		  mol.Molecule(molecule_file='Altrose.pcp', warning_level=1, position=[-5,0,0], scale=400),]


mol = mols[0]
atoms = mol.atoms
# print(mol.torsion_angle(atoms[2], atoms[0], atoms[1], atoms[7]))

# print(id(atoms[0]))

# mm = ff.MM2()
# print(mm.get_energy(mol))

#screen setup
WIDTH, HEIGHT = SIZE = (1200, 720)
screen = scr.Screen3D(SIZE, camera_position=[0., 0, 20.], camera_orientation=(0,0,0), bkgr_colour=(100, 100, 190))
pg.display.set_caption(', '.join([m.name.capitalize() for m in mols]))
clock = pg.time.Clock()
FPS = 120
run = True

#main loop
tick = clock.tick_busy_loop
updt = 0
time = 0

selected_atoms = set()

rot = np.array([0.,0.,0.])
zoom = 0

pg.key.set_repeat()

while run:
	# start = perf_counter()
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT

	screen.draw_axes(100)
	
	screen.clear()
	[screen.draw_shape(m, draw_atoms=True, draw_bonds=True, draw_hydrogens=True) for m in mols]
	[m.rotate(rot) for m in mols]


	keys = pg.key.get_pressed()
	ev = pg.event.get()

	mbd = pg.MOUSEBUTTONDOWN in [e.type for e in ev]

	for e in ev:

		if e.type == pg.MOUSEBUTTONDOWN:
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
			elif e.button == 5:
				if screen.camera_position[2] < 40:
					zoom = dT * screen.camera_position[2] * 3



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

			screen.display_text(f'  {atom.element}{index}  ', (0,0))

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

			screen.display_text(f'  {atoms[0].element}{index[0]}, {atoms[1].element}{index[1]}: {round(atoms[0].distance_to(atoms[1]), 3)} (A)  ', (0,0))

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
				screen.display_text(f'  {atoms[1].element}{index[1]}, {atoms[0].element}{index[0]}, {atoms[2].element}{index[2]}: {round(mols[0].bond_angle(a2, a1, a3, in_degrees=True), 1)} (deg)  ', (0,0))
			elif a1 in a2.bonds and a3 in a2.bonds:
				screen.display_text(f'  {atoms[0].element}{index[0]}, {atoms[1].element}{index[1]}, {atoms[2].element}{index[2]}: {round(mols[0].bond_angle(a1, a2, a3, in_degrees=True), 1)} (deg)  ', (0,0))
			elif a1 in a3.bonds and a2 in a3.bonds:
				screen.display_text(f'  {atoms[0].element}{index[0]}, {atoms[2].element}{index[2]}, {atoms[1].element}{index[1]}: {round(mols[0].bond_angle(a1, a3, a2, in_degrees=True), 1)} (deg)  ', (0,0))
			else:
				screen.display_text(f'  {atoms[0].element}{index[0]}, {atoms[1].element}{index[1]}, {atoms[2].element}{index[2]}  ', (0,0))


	screen.update()

	if keys[pg.K_ESCAPE]:
		run = False
	for event in pg.event.get():
		if event.type == pg.QUIT:
			run = False
			break

	# print(perf_counter()-start)