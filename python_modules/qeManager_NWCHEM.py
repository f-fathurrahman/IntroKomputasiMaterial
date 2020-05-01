from ase import Atoms

def read_nwchem_version(filename):

    f = open(filename, "r")

    while True:
        line = f.readline()
        if not line:
            break
        # Search for string that describes the version
        if "Northwest Computational Chemistry Package (NWChem)" in line:
            version = line.split()[-1]
            break

    f.close()

    return version



def read_nwchem_geoopt(filename):

    f = open(filename, "r")

    while True:
        line = f.readline()
        if not line:
            break

    f.close()

