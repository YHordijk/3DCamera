import math, os, sys, time
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"
import modules.screen4 as scr 
import modules.molecule6 as mol6 
import modules.basisset6 as bs
import modules.colour_maps as cmap
import modules.uff as uff
import modules.utils as utils
import modules.minimizer as minimizer
import modules.plot as plot
import numpy as np
from time import perf_counter
import pubchempy as pcp
import pygame as pg










####### setup
# molecule 				= os.getcwd() + r'\molecules\anthracene.xyz'
molecule 				= 'anthracene'
basis_set 				= 'STO-2G'
repeats 				= 1
add_hydrogens			= False

resolution 				= (1200, 720)
background_colour 		= (100, 100, 190)
draw_axes 				= False
wireframe_mode			= False

do_huckel				= False
pre_render_densities 	= False
points 					= 7000
colour_map 				= cmap.BlueBlackRed(posneg_mode=True)

fancy_format_colours 	= False
fancy_format_time		= True
fancy_format_source		= True
verbosity 				= 0

minimize_structure		= True
sampling_freq			= 1
randomize_structure		= True
plot_energy				= True
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



	utils.message(f'Welcome, this project so far supports the following:', colour='blue')
	[utils.message(f'-- {f}', colour='blue') for f in features]
	utils.message(f'Planned features:', colour='blue')
	[utils.message(f'-- {f}', colour='blue') for f in planned_features]
	print()


MOTD()


utils.ff_print_source(fancy_format_source)
utils.ff_use_colours(fancy_format_colours)
utils.ff_print_time(fancy_format_time)
utils.ff_verbosity(verbosity)

if colour_map == None:
	colour_map = cmap.BlueRed(posneg_mode=True)


#viewer setup
pg.init()
screen = scr.Screen3D(resolution, bkgr_colour=background_colour)
mol = mol6.Molecule(molecule, basis_set_type=basis_set, repeat=repeats)

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
updt_delta = 0
time = 0
rot = np.array([0.,0.])
zoom = 0
pg.key.set_repeat()
run = True
mo_numb = 0
draw_dens = False
camera_range = max(max([a.coords[0] for a in mol.atoms]), max([a.coords[1] for a in mol.atoms]), max([a.coords[2] for a in mol.atoms])) + 10 

screen.camera_position = np.asarray((0,0,camera_range))

utils.message('Please press ENTER to toggle orbital display. Use arrow-keys to switch between orbitals.')
utils.message('Hold CTRL and use mouse to rotate and move molecule.')

#####################
mols = [mol]

if randomize_structure:
	for a1, a2 in mol.get_rotatable_bonds():
		mol.rotate_bond(a1,a2,np.random.random()*2*1*3.14 -1*3.14)

	mol.shake(0.5)
mol.center()
if minimize_structure: 
	mols, energies = minimizer.minimize(mol, 'uff', max_steps=400, sample_freq=5, use_torsions=True)
	if plot_energy:
		p = plot.Plot()
		p.plot(np.arange(len(energies)), energies)
		p.x_label = 'Sample'
		p.y_label = 'Energy (kcal/mol)'
		p.title = f'Progress of Energy Minimization'
		p.show()

screen.size = (resolution)
#####################

drawmol = mols[-1]

while run:
	#tick prep
	updt += 1
	updt_delta += 1
	dT = tick(FPS)/1000
	time += dT

	keys = pg.key.get_pressed()
	ev = pg.event.get()

	screen.clear()


	######################

	if updt_delta == sample_freq*2:
		drawmol = mols[(mols.index(drawmol)+1)%len(mols)]
		updt_delta = 0

	######################


	if draw_dens: screen.draw_density(mos[mo_numb%len(mos)], points, colour_map=colour_map)
	screen.draw_shape(drawmol, wireframe=wireframe_mode, draw_atoms=True, draw_bonds=True, draw_hydrogens=True)	
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
						drawmol.reset_colours()
						for a in drawmol.atoms:
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
					
	[m.rotate(rot) for m in mols]
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
	try:
		l = len(selected_atoms)
		if l > 0:
			sa = list(selected_atoms)
			if l == 1:
				atom = list(selected_atoms)[0]

				if atom in drawmol.atoms:
					index = drawmol.atoms.index(atom)
				else:
					index = drawmol.atoms.index(atom)

				screen.display_text(f'  {atom.symbol}{index+1} {atom.coords} ', (10,10))

			elif l == 2:
				atoms = list(selected_atoms)
				index = []
				for atom in atoms:
					if atom in drawmol.atoms:
						index.append(drawmol.atoms.index(atom))
					else:
						index.append(drawmol.atoms.index(atom))

				screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[1].symbol}{index[1]+1}: {round(atoms[0].distance_to(atoms[1]), 3)} (A)  ', (10,10))

			elif l == 3:
				atoms = list(selected_atoms)
				index = []
				for atom in atoms:
					if atom in drawmol.atoms:
						index.append(drawmol.atoms.index(atom))
					else:
						index.append(drawmol.atoms.index(atom))

				a1, a2, a3 = atoms
				if a2 in a1.bonds and a3 in a1.bonds:
					screen.display_text(f'  {atoms[1].symbol}{index[1]+1}, {atoms[0].symbol}{index[0]+1}, {atoms[2].symbol}{index[2]+1}: {round(drawmol.bond_angle(a2, a1, a3, in_degrees=True), 1)} (deg)  ', (10,10))
				elif a1 in a2.bonds and a3 in a2.bonds:
					screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[1].symbol}{index[1]+1}, {atoms[2].symbol}{index[2]+1}: {round(drawmol.bond_angle(a1, a2, a3, in_degrees=True), 1)} (deg)  ', (10,10))
				elif a1 in a3.bonds and a2 in a3.bonds:
					screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[2].symbol}{index[2]+1}, {atoms[1].symbol}{index[1]+1}: {round(drawmol.bond_angle(a1, a3, a2, in_degrees=True), 1)} (deg)  ', (10,10))
				else:
					screen.display_text(f'  {atoms[0].symbol}{index[0]+1}, {atoms[1].symbol}{index[1]+1}, {atoms[2].symbol}{index[2]+1}  ', (10,10))
	except:
		pass

	screen.update()
