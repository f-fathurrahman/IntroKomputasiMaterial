from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.units import Bohr, Hartree
import numpy as np

atoms = bulk("Si", "diamond", a=10.2631*Bohr)

kden = 2.5

for ec in np.linspace(10.0, 50.0, 9):

    ecutwfc = ec*Hartree

    calc = GPAW( mode=PW(ecutwfc),
                 kpts={ "density": kden },
                 txt="LOG_Si_fcc_ecutwfc_" + str(ec) )
    atoms.set_calculator(calc)

    etot = atoms.get_potential_energy()
    print("%6.3f  %18.10f  %18.10f" % (ec, etot, etot/Hartree))
