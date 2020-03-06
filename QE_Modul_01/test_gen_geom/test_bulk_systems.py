from bulk_systems import *

atoms = build_bulk_fcc("H", 1.5)
atoms.write("GEOM_H_fcc.xsf")

atoms = build_bulk_fcc_cubic("H", 1.5)
atoms.write("GEOM_H_fcc_cubic.xsf")

