import os
import sys

from ase.build import bulk
import numpy as np

from qeManager_PWSCF import *

atoms = bulk("Si")

pspFiles = ['Si.pbe-n-kjpaw_psl.0.1.UPF']

pwinput = PWSCFInput(atoms, pspFiles, filename='PWINPUT',
            kpt_automatic=True, Nk=[8,8,8])

pwinput.filename = 'PWINPUT_scf'
pwinput.CONTROL.pseudo_dir = '/home/efefer/pseudo'
pwinput.set_smearing()
pwinput.write()
os.system('pw.x < PWINPUT_scf > LOG_scf')


pwinput.set_calc_bands('fcc',Nkpts=100)
pwinput.filename = 'PWINPUT_bands'
pwinput.write()
os.system('pw.x < PWINPUT_bands > LOG_bands')

write_bands_inp()
os.system('bands.x < bands.inp > LOG_collect_bands')


Nbands, Nkpts = get_Nbands_Nkpts('LOG_bands')
print('Nbands = ', Nbands, 'Nkpts = ', Nkpts)

#xcoords = pwinput.bands_xcoords
#ebands = read_bandstructure('LOG_bands')

xcoords, ebands = read_bands_gnu('bands.out.gnu', Nbands, Nkpts)

Efermi = read_fermi_level('LOG_scf')

Emin = np.min(ebands)
Emax = np.max(ebands)

print('Emin   = ', Emin)
print('Emax   = ', Emax)
print('Efermi = ', Efermi)

import matplotlib.pyplot as plt
plt.clf()
for ib in range(Nbands):
    plt.plot( xcoords[ib,:], ebands[ib,:], marker='o' )

plt.grid()
plt.ylim(Emin,Emax)

#specialk_xcoords = pwinput.specialk_xcoords
specialk_xcoords = read_specialk_xcoords('LOG_collect_bands')

kmin = np.min(specialk_xcoords)
kmax = np.max(specialk_xcoords)

plt.xlim(kmin,kmax)

for p in specialk_xcoords:
  plt.plot([p, p], [Emin, Emax], 'k-')

plt.plot([kmin,kmax],[Efermi,Efermi], 'k-', linewidth=2)

plt.savefig('bands.png', dpi=300)
