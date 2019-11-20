from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
from ase.units import Bohr, Hartree

atoms = bulk("Si", "diamond", a=10.2631*Bohr)
atoms.write("STRUCT_Si.xsf")

ecutwfc = 15.0*Hartree
calc = GPAW( mode=PW(ecutwfc), 
             setups="hgh",
             xc="LDA_X+LDA_C_VWN",
             eigensolver="dav",
             spinpol=False,
             occupations=FermiDirac(0),
             kpts={"size": (3, 3, 3)},
             txt="-")
atoms.set_calculator(calc)

e1 = atoms.get_potential_energy()
calc.write("TEMP_Si.gpw")

print("===================================================")

print("Total energy      : %18.10f eV = %18.10f Hartree" % (e1,e1/Hartree))

print("===================================================")

