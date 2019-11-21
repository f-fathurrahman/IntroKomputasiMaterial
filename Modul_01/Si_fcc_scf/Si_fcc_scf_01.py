from ase.build import bulk
from gpaw import GPAW, PW
from ase.units import Bohr, Hartree

atoms = bulk("Si", "diamond", a=10.2631*Bohr)

ecutwfc = 15.0*Hartree
calc = GPAW( mode=PW(ecutwfc), 
             kpts={ "size": (3, 3, 3) }, # or simply use kpts=(3,3,3)
             txt="-" )
atoms.set_calculator(calc)

e1 = atoms.get_potential_energy()
calc.write("RESTART_Si_fcc_01.gpw")

print("Total energy : %18.10f eV = %18.10f Hartree" % (e1,e1/Hartree))

