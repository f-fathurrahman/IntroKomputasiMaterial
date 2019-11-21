import numpy as np
import matplotlib.pyplot as plt
import math
import ase.io

def do_plot(prefix):

    print("prefix = " + prefix)

    dEtot = None
    Etot_old = 0.0

    x_kden = []
    y = []

    for i, kden in enumerate( np.linspace(1.0, 5.0, 9) ):  # FIXME: use glob
        filename = prefix + "_" + str(kden)
        res = ase.io.read(filename)
        Etot = res.get_total_energy()

        x_kden.append(kden)
        y.append(Etot)

        if i != 0:
            dEtot = abs(Etot - Etot_old)
            print("%18.10f %18.10f %18.10e" % (kden, Etot, dEtot))
        else:
            print("%18.10f %18.10f" % (kden, Etot))
        Etot_old = Etot

    y = np.array(y)
    y = y - y[-1]
    plt.clf()
    plt.plot(x_kden, y, marker="o")
    plt.grid()
    plt.savefig("IMG_" + prefix + ".pdf")


do_plot("LOG_kden")
do_plot("LOG_kden_v2")
do_plot("LOG_Cu_fcc_kden")