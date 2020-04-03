"""

"""

import sys
sys.path.append("..")

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
import matplotlib.cm as mc
import matplotlib.colors as mcol

import matplotlib as mpl

from analysis.plotting import *
import numpy as np


def forcing(z):
    return -np.cos(z)

def c_inh(z):
    return -1+ 1/(1.5+forcing(z))

def main():
    z = np.linspace(-np.pi,  np.pi , 100)

    n = 5
    arctan_args = np.linspace(0, n-1, n*3)

    norm = mpl.colors.Normalize(vmin=arctan_args.min(),
                                vmax=arctan_args.max())

    cmap = mpl.cm.ScalarMappable(norm=norm, cmap='copper')
    cmap.set_array([])

    fig, ax = plt.subplots()

    for arg_i in arctan_args:
        delta_phi = np.arctan( arg_i )
        ax.plot(z, np.cos(z+delta_phi), c=cmap.to_rgba(arg_i))

    ax.plot(z, forcing(z), label = r'contractile stress' , c='darkblue')
    # ax.plot(z, c_inh(z), label = r'$c = \frac{1}{const. + f_i} -1 $' , c='darkgreen')

    ax.set_xlabel(r'$z = x/\lambda  - t \nu$')
    # ax2.legend()
    ticks = np.array([-np.pi, -np.pi/2, 0, np.pi/2 ,np.pi])
    labels = [r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    cbar = plt.colorbar(cmap)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(r'radius for $2 \kappa + \gamma$', rotation=270)

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center', ncol=n+1)

    plt.show()


if __name__ == "__main__":
    main()
