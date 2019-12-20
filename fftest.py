import modules.molecule3 as mol 
import modules.forcefields as forcefields
import os

mols = [mol.Molecule(molecule_file='ethane.pcp', warning_level=0),
		mol.Molecule(molecule_file=f'Molecules\\ethane eclipsed.xyz', warning_level=0),]
		# mol.Molecule(molecule_file='ethane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='propane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='butane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='pentane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='hexane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='heptane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='octane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='isobutane.pcp', warning_level=0),
		# mol.Molecule(molecule_file='isopentane.pcp', warning_level=0)]


# ff = forcefields.MM2()

# energies = [ff.get_energy(mol) for mol in mols]

# print(energies)

for t in mols[0].get_unique_torsion_angles(in_degrees=True):
	print(t)