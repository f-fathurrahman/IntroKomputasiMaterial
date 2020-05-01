from qeManager_NWCHEM import *

version = read_nwchem_out_version("../WORKS_NWCHEM/C2H2O/MOLECULE.nwout")
print("version = ", version)