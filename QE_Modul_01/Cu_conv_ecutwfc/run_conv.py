from ase.build import bulk
from ase.units import Bohr, Rydberg
import ase.io

import sys
sys.path.append("../")
from qeManager_PWSCF import *

atoms = bulk("Cu", "fcc")

pspFiles = ["Cu.pbe-kjpaw.UPF"]

pwinput = PWSCFInput(atoms, pspFiles, filename="PWINPUT", kpt_automatic=True, Nk=[8,8,8])
pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo"
pwinput.SYSTEM.occupations = "smearing"
pwinput.SYSTEM.smearing = "mv"
pwinput.SYSTEM.degauss = 0.01
conv_test = ConvergenceTest( pwinput, what="ecutwfc", values=np.linspace(20.0, 100.0, 9) )

#conv_test.run()

ecutrho, energies, forces = conv_test.read()

Ndata = len(ecutrho)
Eref = energies[-1]  # take the most converged data (the last one) as the reference

x_ecutwfc = []
y = []
print("\nConvergence in total energy: (eV)")
for i in range(Ndata):
    ec = ecutrho[i]
    ene = energies[i]*Rydberg
    dEtot = (energies[i]-Eref)*Rydberg
    x_ecutwfc.append(ec)
    y.append(ene)
    print("%10.5f %18.10f %18.10f" % (ec, ene, dEtot))

print("\nConvergence in total force:")
for i in range(Ndata):
    totForces = np.sum(np.abs(forces[i]))/len(atoms)
    print("%10.5f %18.10f" % (ecutrho[i], totForces))


prefix = "Cu_fcc_ecutwfc"

y = np.array(y)
y = y - y[-1]

import matplotlib.pyplot as plt
plt.clf()
plt.plot(x_ecutwfc, y, marker="o")
plt.ylim(-0.001, 0.05)
plt.ylabel("Relative energy [eV]")
plt.xlabel("ecutwfc [Ry]")
plt.grid()
plt.savefig("IMG_" + prefix + ".pdf")
