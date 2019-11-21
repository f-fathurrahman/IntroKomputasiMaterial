import numpy as np
import matplotlib.pyplot as plt
import math
import ase.io

def do_plot(prefix):

    print("prefix = " + prefix)

    dEtot = None
    Etot_old = 0.0

    x_ec = []
    y = []

    for i, ec in enumerate( np.linspace(10.0, 50.0, 9) ):  # FIXME: use glob
        filename = prefix + "_" + str(ec)
        res = ase.io.read(filename)
        Etot = res.get_total_energy()

        x_ec.append(ec)
        y.append(Etot)

        if i != 0:
            dEtot = -(Etot - Etot_old)
            print("%18.10f %18.10f %18.10e" % (ec, Etot, dEtot))
        else:
            print("%18.10f %18.10f" % (ec, Etot))
        Etot_old = Etot

    y = np.array(y)
    y = y - y[-1]
    plt.clf()
    plt.plot(x_ec, y, marker="o")
    plt.ylim(-0.001, 0.05)
    plt.grid()
    plt.savefig("IMG_" + prefix + ".pdf")


do_plot("LOG_Cu_fcc_ecutwfc")