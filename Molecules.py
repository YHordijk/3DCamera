import math, os, sys, time
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"
import modules.screen4 as scr 
import modules.molecule6 as mol6 
import modules.basisset6 as bs
import modules.colour_maps as cmap
import modules.uff as uff
import modules.utils as utils
import numpy as np
from time import perf_counter
import pubchempy as pcp
import pygame as pg









####### setup
molecule 				= 'cyclohexane'
basis_set 				= 'STO-2G'
repeats 				= 1
add_hydrogens			= False

resolution 				= (1200, 720)
background_colour 		= (100, 100, 190)
draw_axes 				= False
wireframe_mode			= True

do_huckel				= False
pre_render_densities 	= False
points 					= 7000
colour_map 				= cmap.BlueBlackRed(posneg_mode=True)

fancy_format_colours 	= False
fancy_format_time		= True
fancy_format_source		= True

calculate_uff_energy 	= True
uff_verbosity			= 0
#######











def MOTD():
	utils.ff_print_source(False)
	utils.ff_use_colours(fancy_format_colours)
	utils.ff_print_time(False)
	os.system('color 07')
	features = ['Loading molecules from xyz files or alternatively download from pubchem.',
				'Visualising molecules as a ball-and-stick model or as a wireframe model.',
				'Algorithms for guessing atom types, bonds, and bond orders.',
				'Support for STO-nG basis set for atomic/molecular orbital visualtions and calculations.',
				'Crude extended-Hückel method for calculation of molecular orbitals.']
	planned_features = ['Improvements to extended-Hückel and molecular integrals.',
						'Marching cubes algorithm to visualise orbitals as iso-volumes instead of current random-dot style.']



	utils.message(f'Welcome, this project so far supports the following:', 'blue')
	[utils.message(f'-- {f}', 'blue') for f in features]
	utils.message(f'Planned features:', 'blue')
	[utils.message(f'-- {f}', 'blue') for f in planned_features]
	print()


MOTD()


utils.ff_print_source(fancy_format_source)
utils.ff_use_colours(fancy_format_colours)
utils.ff_print_time(fancy_format_time)

if colour_map == None:
	colour_map = cmap.BlueRed(posneg_mode=True)


#viewer setup
pg.init()
screen = scr.Screen3D(resolution, bkgr_colour=background_colour)
mol = mol6.Molecule(molecule, basis_set_type=basis_set, repeat=repeats)

if calculate_uff_energy:
	ff = uff.ForceField()
	ff.get_energy(mol, verbosity=uff_verbosity)

if add_hydrogens:
	mol.add_hydrogens()

atoms = mol.atoms
pg.display.set_caption(mol.name)
selected_atoms = set()




bs.extended_huckel(mol) 
mos = mol.molecular_orbitals

if pre_render_densities: screen.pre_render_densities(mos, points=points, colour_map=colour_map)

clock = pg.time.Clock()
FPS = 120
tick = clock.tick_busy_loop
updt = 0
time = 0
rot = np.array([0.,0.])
zoom = 0
pg.key.set_repeat()
run = True
mo_numb = 0
draw_dens = False
camera_range = max(max([a.coords[0] for a in mol.atoms]), max([a.coords[1] for a in mol.atoms]), max([a.coords[2] for a in mol.atoms])) + 3 

screen.camera_position = np.asarray((0,0,camera_range))

utils.message('Please press ENTER to toggle orbital display. Use arrow-keys to switch between orbitals.')
utils.message('Hold CTRL and use mouse to rotate and move molecule.')

##############


b = []
for bond in mol.get_unique_bonds():
	b.append([*bond, bond[0].distance_to(bond[1]), np.random.normal(0,math.pi)])

angles = []
for a1, a2, a3, a4, phi in mol.get_unique_torsion_angles():
	angles.append([a1,a2,a3,a4,phi,np.random.normal(-math.pi,math.pi)/1000])

##############


while run:
	#tick prep
	updt += 1
	dT = tick(FPS)/1000
	time += dT

	keys = pg.key.get_pressed()
	ev = pg.event.get()

	screen.clear()


	################ 


	for a1, a2, r, o in b:
		mol.stretch_bond(a1, a2, r + math.sin(10*time + o)/10)
	for a1, a2, a3, a4, phi, o in angles:
		# mol.rotate_bond(a2, a3, o/2)
		mol.rotate_bond(a3, a2, o*2)

	mol.center()
	# mol.align_bond_to_vector(mol.atoms[1], mol.atoms[0], (math.sin(time),math.cos(time),0))

	################


	if draw_dens: screen.draw_density(mos[mo_numb%len(mos)], points, colour_map=colour_map)
	screen.draw_shape(mol, wireframe=wireframe_mode, draw_atoms=True, draw_bonds=True, draw_hydrogens=True)	
	if draw_axes: screen.draw_axes(1)

	if keys[pg.K_ESCAPE]:
		run = False

	for e in ev:
		if e.type == pg.VIDEORESIZE:
			screen.size = e.dict['size']

		if e.type == pg.QUIT:
			run = False

		if e.type == pg.KEYDOWN:
			if e.key == pg.K_RIGHT:
				mo_numb += 1
			if e.key == pg.K_LEFT:
				mo_numb -= 1
			if e.key == pg.K_RETURN:
				draw_dens = 1 - draw_dens

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
					
	mol.rotate(rot)
	rot *= 0.8
	move = pg.mouse.get_rel()
	if keys[pg.K_LCTRL] or keys[pg.K_RCTRL]:
		if pg.mouse.get_pressed()[2]:
				screen.camera_position[0] += move[0]/50
				screen.camera_position[1] += move[1]/50

		if pg.mouse.get_pressed()[0]:
			rot = np.asarray([move[1]/150, -move[0]/150])


	screen.camera_position[2] += zoom
	zoom *= 0.8


	#tick end
	l = len(selected_atoms)
	if l > 0:
		sa = list(selected_atoms)
		if l == 1:
			atom = list(selected_atoms)[0]

			if atom in mol.atoms:
				index = mol.atoms.index(atom)
			else:
				index = mol.atoms.index(atom)

			screen.display_text(f'  {atom.symbol}{index+1} {atom.coords} ', (10,10))

		elif l == 2:
			atoms = list(selected_atoms)
			index = []
			for atom in atoms:
				if atom in mol.atoms:
					index.append(mol.atoms.index(atom))
				else:
					index.append(mol.atoms.index(atom))

			screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[1].symbol}{index[1]+1}: {round(atoms[0].distance_to(atoms[1]), 3)} (A)  ', (10,10))

		elif l == 3:
			atoms = list(selected_atoms)
			index = []
			for atom in atoms:
				if atom in mol.atoms:
					index.append(mol.atoms.index(atom))
				else:
					index.append(mol.atoms.index(atom))

			a1, a2, a3 = atoms
			if a2 in a1.bonds and a3 in a1.bonds:
				screen.display_text(f'  {atoms[1].symbol}{index[1]+1}, {atoms[0].symbol}{index[0]+1}, {atoms[2].symbol}{index[2]+1}: {round(mol.bond_angle(a2, a1, a3, in_degrees=True), 1)} (deg)  ', (10,10))
			elif a1 in a2.bonds and a3 in a2.bonds:
				screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[1].symbol}{index[1]+1}, {atoms[2].symbol}{index[2]+1}: {round(mol.bond_angle(a1, a2, a3, in_degrees=True), 1)} (deg)  ', (10,10))
			elif a1 in a3.bonds and a2 in a3.bonds:
				screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[2].symbol}{index[2]+1}, {atoms[1].symbol}{index[1]+1}: {round(mol.bond_angle(a1, a3, a2, in_degrees=True), 1)} (deg)  ', (10,10))
			else:
				screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[1].symbol}{index[1]+1}, {atoms[2].symbol}{index[2]+1}  ', (10,10))


	screen.update()
