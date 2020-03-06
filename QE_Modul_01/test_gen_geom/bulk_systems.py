# Various functions to build common bulk systems
# Stripped down from ASE

from ase import Atoms

def build_bulk_fcc(name, a):
    b = a / 2
    atoms = Atoms(name, cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)
    return atoms

def build_bulk_fcc_cubic(name, a):
    atoms = Atoms(4*name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0.5, 0.5, 0)])
    return atoms

def build_bulk_bcc(name, a):
    b = a / 2
    atoms = Atoms(name, cell=[(-b, b, b), (b, -b, b), (b, b, -b)], pbc=True)
    return atoms


def build_bulk_bcc_cubic(name, a):
    atoms = Atoms(2*name, cell=(a, a, a), pbc=True,
                        scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    return atoms
