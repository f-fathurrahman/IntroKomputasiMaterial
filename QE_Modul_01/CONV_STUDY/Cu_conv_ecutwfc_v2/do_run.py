import sys
sys.path.append("/home/efefer/WORKS/my_github_repos/IntroKomputasiMaterial/python_modules")

from qeManager_PWSCF import *

from ase.build import bulk
from ase.units import Bohr, Rydberg
import ase.io

MODE = 2 #
# 0: generate input file
# 1: run
# 2:

atoms = bulk("Cu", "fcc")

pspFiles = ["Cu_ONCV_PBE-1.0.upf"]

ECUTLIST = np.arange(50.0, 220.0, 10.0) # cutoff energy in Ry
Ndata = len(ECUTLIST)

# Generate pwinput objects
INP_FILELIST = []
OUT_FILELIST = []
PWINPUTS = []
for i, ec in enumerate(ECUTLIST):
    #
    pwinput = PWSCFInput(atoms, pspFiles, filename="PWINPUT", kpt_automatic=True, Nk=[8,8,8])
    pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo/ONCV_PBE"
    pwinput.SYSTEM.ecutwfc = ec
    pwinput.SYSTEM.ecutrho = 4.0*ec
    pwinput.set_smearing()
    #
    pwinput.filename = "INP_PW_ecutwfc_" + str(ec)
    INP_FILELIST.append(pwinput.filename)
    PWINPUTS.append(pwinput)
    #
    fileout = "OUT_PW_ecutwfc_" + str(ec)
    OUT_FILELIST.append(fileout)

# Write input files
if MODE == 0:
    for pwinput in PWINPUTS:
        pwinput.write()

if MODE == 1:
    for filein, fileout in zip(INP_FILELIST, OUT_FILELIST):
        os.system("pw.x < " + filein + " > " + fileout)
        print("Running " + filein)

energies = []
if MODE == 2:

    for fileout in OUT_FILELIST:
        energy = read_pwscf_energy(fileout)
        energies.append(energy)
    
    Eref = energies[-1] # - 1e-12 # Add small eps to avoid log(0)
    y = []
    print("------------------------------------------------")
    print("ecutwfc (Ry)      Etot (eV)           ΔEtot (eV)")
    print("------------------------------------------------")
    for i in range(Ndata):
        ec = ECUTLIST[i]
        ene = energies[i]
        dEtot = ( energies[i] - Eref )
        print("%10.5f %18.10f %18.10f" % (ec, ene*Rydberg, dEtot*Rydberg))
        y.append(dEtot*Rydberg)

    import matplotlib.pyplot as plt
    import numpy as np

    plt.clf()
    plt.plot(ECUTLIST[:Ndata-1], y[:Ndata-1], marker="o")
    plt.grid()
    plt.xlabel("ecutwfc (Ry)")
    plt.ylabel("ΔEtot (eV)")
    plt.ylim(0.0, 0.05)
    plt.title(pspFiles[0])
    plt.savefig("IMG_conv_ecutwfc_lin.png", dpi=150)
    plt.savefig("IMG_conv_ecutwfc_lin.pdf")

    plt.clf()
    plt.plot(ECUTLIST[:Ndata-1], np.log10(y[:Ndata-1]), marker="o")
    plt.grid()
    plt.xlabel("ecutwfc (Ry)")
    plt.ylabel("log10(ΔEtot (eV))")
    plt.title(pspFiles[0])
    plt.savefig("IMG_conv_ecutwfc_log.png", dpi=150)    
    plt.savefig("IMG_conv_ecutwfc_log.pdf")
