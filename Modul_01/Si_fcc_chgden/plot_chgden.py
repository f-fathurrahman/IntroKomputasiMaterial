from gpaw import GPAW, restart
import matplotlib.pyplot as plt
import ase.io

atoms, calc = restart("../Si_fcc_scf/RESTART_Si_fcc_01.gpw")
gridrefinement = 2.0  # this is the maximum value that is implemented currently

# Get pseudo and all-electron densities
rhoe_PS = calc.get_pseudo_density()
rhoe_AE = calc.get_all_electron_density(gridrefinement=gridrefinement)

import numpy as np
np.save("rhoe_PS.dat", rhoe_PS)
np.save("rhoe_AE.dat", rhoe_AE)

from xsf_write import *

xsf_write("Si_rhoe_PS.xsf", atoms, data=rhoe_PS)
xsf_write("Si_rhoe_AE.xsf", atoms, data=rhoe_AE)

# This is buggy (at the moment)
#ase.io.write("Si_rhoe_PS.xsf", atoms, format="xsf", data=rhoe_PS_p)
#ase.io.write("Si_rhoe_AE.xsf", atoms, format="xsf", data=rhoe_AE)