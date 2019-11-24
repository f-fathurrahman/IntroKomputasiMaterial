import sys
sys.path.append("../../python_modules")

from qeManager_PWSCF import *

from ase.build import molecule
from ase.units import Bohr, Rydberg
import ase.io

atoms = molecule("N2")
atoms.set_pbc([True,True,True])
atoms.center(vacuum=5.0)

atoms.write("N2.xsf")

pspFiles = ["N.pbe-n-kjpaw_psl.0.1.UPF "]

pwinput = PWSCFInput(atoms, pspFiles, filename="PWINPUT", gamma_only=True)
pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo/PSLIB"

ECUTLIST = np.linspace(10.0, 100.0, 10) # cutoff energy in Ry

DUAL = 4

# write files
INP_FILELIST = []
for i, ec in enumerate(ECUTLIST):
    
    pwinput.SYSTEM.ecutwfc = ec
    pwinput.SYSTEM.ecutrho = DUAL*ec
    pwinput.SYSTEM.assume_isolated = "mt"

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


ene = np.zeros(len(OUT_FILELIST))
for i, f in enumerate(OUT_FILELIST):
    ene[i] = read_pwscf_energy(f)
    print("ene[%d] = %18.10f" % (i, ene[i]))

Eref = ene[-1]  # take the most converged data (the last one) as the reference
ene[:] = (ene[:] - Eref)*Ry # convert to eV

import matplotlib.pyplot as plt
plt.clf()
plt.plot(ECUTLIST, ene, marker="o")
plt.ylim(-0.001, 0.05)
plt.ylabel("Relative energy [eV]")
plt.xlabel("ecutwfc [Ry]")
plt.grid()
plt.savefig("IMG_conv_ecutwfc.pdf")
