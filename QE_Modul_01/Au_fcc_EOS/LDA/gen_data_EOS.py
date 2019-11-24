import sys
sys.path.append("../../../python_modules")

from qeManager_PWSCF import *

from ase.build import bulk
import ase.io
from ase.units import Bohr, kJ, Ry
from ase.eos import EquationOfState

pspFiles = ["au_lda_v1.uspp.F.UPF"]

a_exp = 4.065
a_min = 0.8*a_exp
a_max = 1.2*a_exp
A_LIST = np.linspace(a_min, a_max, 10)

# write files
INP_FILELIST = []
for i, a in enumerate(A_LIST):
    
    atoms = bulk("Au", "fcc", a=a)
    
    pwinput = PWSCFInput(atoms, pspFiles, filename="PWINPUT", kpt_automatic=True, Nk=[10,10,10])

    pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo/GBRV_LDA"
    pwinput.SYSTEM.ecutwfc = 40.0
    pwinput.SYSTEM.ecutrho = 200.0

    pwinput.SYSTEM.occupations = "smearing"
    pwinput.SYSTEM.smearing = "mv"
    pwinput.SYSTEM.degauss = 0.1/Ry  # 0.1 eV to Ry
    
    pwinput.filename = "INP_PW_" + str(i)
    
    INP_FILELIST.append(pwinput.filename)
    pwinput.write()

# run commands
OUT_FILELIST = []
for i, f in enumerate(INP_FILELIST):
    fileout = "OUT_PW_" + str(i) 
    OUT_FILELIST.append(fileout)
    # Uncomment this if the run is already done
    os.system("pw.x < " + f + " > " + fileout)
    print("Running " + str(i) + " is done")


ene = np.zeros(len(OUT_FILELIST))
for i, f in enumerate(OUT_FILELIST):
    ene[i] = read_pwscf_energy(f)
    print("ene[%d] = %18.10f" % (i, ene[i]))


volumes = A_LIST[:]**3/4.0 # fcc volumes (in Angstrom^3)
ene[:] = ene[:]*Ry # convert to eV

eos = EquationOfState(volumes, ene)
v0, e0, B = eos.fit()
print("B = %18.10f GPa" % (B/kJ * 1.0e24))
eos.plot("IMG_fit_EOS.pdf")

a_calc = (4*v0)**(1/3)
print("a_exp  = ", a_exp)
print("a_calc = ", a_calc)