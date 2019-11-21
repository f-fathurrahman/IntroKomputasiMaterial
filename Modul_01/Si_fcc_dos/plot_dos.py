from gpaw import GPAW, restart
import matplotlib.pyplot as plt

atoms, calc = restart('TEMP_Si_k8.gpw')

# Restart from ground state and fix potential:
calc = GPAW('Si_gs.gpw',
            nbands=16,
            fixdensity=True,
            symmetry='off',
            kpts={'path': 'GXWKL', 'npoints': 60},
            convergence={'bands': 8})

calc.get_potential_energy()

x_energy, dos = calc.get_dos(spin=0, npts=2000, width=0.2)
E_fermi = calc.get_fermi_level()

plt.clf()
plt.plot(x_energy - E_fermi, dos)
#plt.axis([-15, 10, None, 4])
plt.ylabel("DOS")
plt.savefig("Si_total_DOS_k8.pdf")

print("Npts = ", len(x_energy))