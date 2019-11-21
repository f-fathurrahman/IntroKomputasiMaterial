from ase.build import bulk

atoms = bulk("Li", "bcc")
atoms.write("STRUCT_Li_bcc.xsf")