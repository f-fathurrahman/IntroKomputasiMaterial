from ase import Atoms

def read_nwchem_out_version(filename):
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

def read_nwchem_out_Natoms(filename):
    Natoms = 0
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        # Search for string that describes the version
        if "XYZ format geometry" in line:
            line = f.readline() # skip
            Natoms = int( f.readline() )
            break
    f.close()
    return Natoms

def read_nwchem_out_geoopt(filename):

    f = open(filename, "r")

    while True:
        line = f.readline()
        if not line:
            break
        if "Geometry " in line:
            for i in range(6):
                line = f.readline()
            line = f.readline()
            print(line)
            break

    f.close()

