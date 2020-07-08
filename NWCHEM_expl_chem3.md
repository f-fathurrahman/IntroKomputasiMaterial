# Formaldehyde energy

Gaussian09
```
Int=(UltraFine,Acc2E=12)
SCF Done:  E(RB3LYP) =  -114.544267470     A.U. after    9 cycles
 Dipole moment (field-independent basis, Debye):
    X=              0.0000    Y=              0.0000    Z=             -2.4337  Tot=              2.4337
```

NWCHEM (using spherical basis)
```
Total DFT energy =     -114.544267468654
   Dipole moment        2.4337307197 Debye(s)
             DMX       -0.0000000000 DMXEFC        0.0000000000
             DMY        0.0000000000 DMYEFC        0.0000000000
             DMZ       -2.4337307197 DMZEFC        0.0000000000
   -EFC- dipole         0.0000000000 DEBYE(S)
   Total dipole         2.4337307197 DEBYE(S)
```



# Comparing formaldehyde and acetone

Result from Gaussian09
```
Using default
SCF Done:  E(RB3LYP) =  -193.215648547     A.U. after   11 cycles
 Dipole moment (field-independent basis, Debye):
X=              0.0000    Y=              0.0000    Z=             -3.2420  Tot=              3.2420

Using Int=(UltraFine,Acc2E=12)
 Dipole moment (field-independent basis, Debye):
SCF Done:  E(RB3LYP) =  -193.215654442     A.U. after   11 cycles
X=              0.0000    Y=              0.0000    Z=             -3.2420  Tot=              3.2420
```

Result from NWCHEM default (spherical basis:
```
Using default
Total DFT energy =     -193.215654438824
   Dipole moment        3.2420084083 Debye(s)
             DMX       -0.0000000000 DMXEFC        0.0000000000
             DMY        0.0000000000 DMYEFC        0.0000000000
             DMZ       -3.2420084083 DMZEFC        0.0000000000
   -EFC- dipole         0.0000000000 DEBYE(S)
   Total dipole         3.2420084083 DEBYE(S)

Using grid xfine; convergence density 1e-8
Total DFT energy =     -193.215654566424
   Dipole moment        3.2420229534 Debye(s)
             DMX       -0.0000000000 DMXEFC        0.0000000000
             DMY        0.0000000000 DMYEFC        0.0000000000
             DMZ       -3.2420229534 DMZEFC        0.0000000000
   -EFC- dipole         0.0000000000 DEBYE(S)
   Total dipole         3.2420229534 DEBYE(S)
```