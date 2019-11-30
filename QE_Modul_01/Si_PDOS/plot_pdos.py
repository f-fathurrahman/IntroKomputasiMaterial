import numpy as np
import matplotlib.pyplot as plt

DIRNAME = "DATA_PDOS_v2"

E_Fermi = 6.3365

data = np.loadtxt(DIRNAME + "/pwscf.pdos_atm#1(Si)_wfc#1(s)")

ene = data[:,0] - E_Fermi

ldos_s = data[:,1]
pdos_s = data[:,2]

data = np.loadtxt(DIRNAME + "/pwscf.pdos_atm#1(Si)_wfc#2(p)")

ldos_p = data[:,1]
pdos_pz = data[:,2]
pdos_px = data[:,3]
pdos_py = data[:,4]

plt.clf()
plt.plot(ene, pdos_s, label="s")
plt.plot(ene, pdos_pz, label="pz")
plt.grid()
plt.xlim(-5.0, 5.0)
plt.savefig("IMG_Si_1_dos.pdf")
