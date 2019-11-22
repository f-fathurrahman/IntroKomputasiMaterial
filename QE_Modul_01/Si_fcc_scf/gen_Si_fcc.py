from ase.build import bulk

atoms = bulk("C", "diamond")
atoms.write("STRUCT_C_diamond.xsf")