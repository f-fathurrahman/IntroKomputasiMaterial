from bulk_systems import *

atoms = build_bulk_fcc("H", 1.5)
atoms.write("GEOM_H_fcc.xsf")

atoms = build_bulk_fcc_cubic("H", 1.5)
atoms.write("GEOM_H_fcc_cubic.xsf")

atoms = build_bulk_bcc("H", 1.5)
atoms.write("GEOM_H_bcc.xsf")

atoms = build_bulk_bcc_cubic("H", 1.5)
atoms.write("GEOM_H_bcc_cubic.xsf")
