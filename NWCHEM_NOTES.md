# General notes

Conventions:

- Extension for NWChem input file: `.nw`
- Extension for NWChem log file: `.nwo` or `.nwout`.

Running NWChem:
```
nwchem MOLECULE.nw
```

Running NWChem (alternative):
```
nwchem MOLECULE.nw > MOLECULE.nw
```

## About `echo`

If `echo` is specified in the input file, then the input file will be echoed
to stdout. This is useful for debugging purpose and/or reproducibility.

## Hartree-Fock 

```
scf
  singlet
  maxiter 128
  noprint "final vector analysis"
end
```


## DFT options

```
DFT
  direct
  grid nodisk
  noprint "final vectors analysis"
  maxiter 200
  xc b3lyp
end
```

The double quotes sign is important for `noprint` directive to have effect. Without
double quotes, NWChem will still print final vectors analysis (symmetries and
coefficients).

Using the option `grid nodisk` can reduce disk usage, however, it results in slower
calculation (if reading disk is not a bottleneck in the calculation).

Getting total energy (in single-point calculation):
```
grep "DFT energy" LOG
```

### Open shell DFT

Use keyword `odft`.

### Convergence tricks for metals

Using `CONVERGENCE rabuck` or `SMEAR 0.01`.

## Calculating properties

```
task property
  dipole
end

task dft property
```

It seems that running `task dft energy property` is the same as `task dft property`.

## Molden

```
property
  vectors PREFIX.movecs
  moldenfile
  molden_norm janpa
end
```

The MOs are read from `PREFIX.movecs`.

## Basis spherical

It seems that specifying the basis to be spherical is essential
to obtain similar energy result to Gaussian09.
