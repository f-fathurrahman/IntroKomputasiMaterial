from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.units import Bohr, Hartree
import numpy as np

atoms = bulk("Si", "diamond", a=10.2631*Bohr)

ecutwfc = 15.0*Hartree

for kden in np.linspace(1.0, 5.0, 9):
    calc = GPAW( mode=PW(ecutwfc),
                 occupations=FermiDirac(0.0),
                 kpts={ "density": kden },
                 txt="LOG_kden_v2_" + str(kden) )
    atoms.set_calculator(calc)

    etot = atoms.get_potential_energy()
    print("%6.3f  %18.10f  %18.10f" % (kden, etot, etot/Hartree))
