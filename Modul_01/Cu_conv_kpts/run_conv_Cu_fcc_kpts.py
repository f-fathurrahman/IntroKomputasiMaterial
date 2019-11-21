from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.units import Bohr, Hartree
import numpy as np

atoms = bulk("Cu", "fcc")

ecutwfc = 30.0*Hartree

for kden in np.linspace(1.0, 5.0, 9):
    calc = GPAW( mode=PW(ecutwfc),
                 kpts={ "density": kden },
                 txt="LOG_Cu_fcc_kden_" + str(kden) )
    atoms.set_calculator(calc)

    etot = atoms.get_potential_energy()
    print("%6.3f  %18.10f  %18.10f" % (kden, etot, etot/Hartree))
