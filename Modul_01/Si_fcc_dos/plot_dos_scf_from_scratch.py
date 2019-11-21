from ase.build import bulk
from gpaw import GPAW, PW
from ase.units import Bohr, Hartree
import matplotlib.pyplot as plt

atoms = bulk("Si", "diamond", a=10.2631*Bohr)

ecutwfc = 15.0*Hartree
calc = GPAW( mode=PW(ecutwfc), 
             kpts=(8, 8, 8),
             txt="-" )
atoms.set_calculator(calc)

e1 = atoms.get_potential_energy()
calc.write("RESTART_Si_fcc_k8x8x8.gpw")

print("Total energy : %18.10f eV = %18.10f Hartree" % (e1,e1/Hartree))

x_energy, dos = calc.get_dos(spin=0, npts=2000, width=0.2)
E_fermi = calc.get_fermi_level()

plt.clf()
plt.plot(x_energy - E_fermi, dos)
plt.ylabel("DOS")
plt.savefig("IMG_Si_total_DOS_scf_from_scratch_k8x8x8.pdf")
