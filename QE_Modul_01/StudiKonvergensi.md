# Studi Konvergensi

Fadjar Fathurrahman

# Pendahuluan

Pada tutorial ini akan dibahas mengenai studi konvergensi.

Dalam melakukan perhitungan struktur elektronik, seringkali kita harus menentukan beberapa
parameter atau input variabel pada kalkulasi kita. Parameter tersebut biasanya terkait
dengan metode numerik atau diskritisasi. Pada perhitungan dengan paket program berbasis
gelombang bidang biasanya ada dua parameter yang harus dipelajari konvergensinya:

- energi cutoff: `ecutwfc` dan `ecutrho` untuk Kohn-Sham orbital dan
  kerapatan elektron. Untuk kasus pseudopotensial norm-conserving kita hanya perlu
  melakukan ini untuk `ecutwfc`.

- kpoint grid: untuk diskritisasi Brillouin zone

# Konvergensi Total energy vs ecutwfc pada Si fcc

Berikut ini adalah script Python sederhana yang dapat digunakan untuk melakukan studi
konvergensi terhadap.

```python
import sys
sys.path.append("../../python_modules")
from qeManager_PWSCF import *
```

```python
from ase.build import bulk
from ase.units import Bohr, Rydberg
import ase.io
```

atoms = bulk("Si", "diamond", a=5.431)

pspFiles = ["Si.pz-hgh.UPF"]

pwinput = PWSCFInput(atoms, pspFiles, filename="PWINPUT", kpt_automatic=True, Nk=[4,4,4])
pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo/HGH"

ECUTLIST = np.linspace(10.0, 100.0, 10) # cutoff energy in Ry

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


