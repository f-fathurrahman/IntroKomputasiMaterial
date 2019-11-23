import sys
sys.path.append("../../python_modules")

from qeManager_PWSCF import *

from ase.build import bulk
from ase.units import Bohr, Rydberg
import ase.io

atoms = bulk("Si", "diamond", a=5.431)

pspFiles = ["Si.pz-hgh.UPF"]

pwinput = PWSCFInput(atoms, pspFiles, filename="PWINPUT", kpt_automatic=True, Nk=[4,4,4])
pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo/HGH"

ECUTLIST = np.linspace(10.0, 100.0, 10) # cutoff energy in Ry

print(ECUTLIST)

# write files
INP_FILELIST = []
for i, ec in enumerate(ECUTLIST):
    pwinput.SYSTEM.ecutwfc = ec
    pwinput.filename = "INP_PW_ecutwfc_" + str(i)
    INP_FILELIST.append(pwinput.filename)
    pwinput.write()

# run commands
OUT_FILELIST = []
for i, f in enumerate(INP_FILELIST):
    fileout = "OUT_PW_ecutwfc_" + str(i) 
    OUT_FILELIST.append(fileout)
    os.system("pw.x < " + f + " > " + fileout)
    print("Running " + str(i) + " is done")


"""
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
"""