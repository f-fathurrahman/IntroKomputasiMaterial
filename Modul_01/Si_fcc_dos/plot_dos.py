from gpaw import GPAW, restart
import matplotlib.pyplot as plt

# Restart from ground state and fix potential:
calc = GPAW( "../Si_fcc_scf/RESTART_Si_fcc_01.gpw",
             fixdensity=True,
             symmetry="off",
             kpts=(8,8,8) )

calc.get_potential_energy()

x_energy, dos = calc.get_dos(spin=0, npts=2000, width=0.2)
E_fermi = calc.get_fermi_level()

plt.clf()
plt.plot(x_energy - E_fermi, dos)
#plt.axis([-15, 10, None, 4])
plt.ylabel("DOS")
plt.savefig("IMG_Si_total_DOS_nscf_k8.pdf")
