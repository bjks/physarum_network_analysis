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


def stress(c, sens, u=1):
    return 3 - sens * c/(c + u)


def main():
    c = np.linspace(0,  2 , 100)

    sensitivity = np.linspace(0.5, 5, 10)

    norm = mpl.colors.Normalize(vmin=sensitivity.min()-3,
                                vmax=sensitivity.max())

    cmap = mpl.cm.ScalarMappable(norm=norm, cmap='Wistia')
    cmap.set_array([])

    fig, ax = plt.subplots()

    for sens_i in sensitivity:
        ax.plot(c, stress(c, sens_i), c=cmap.to_rgba(sens_i))

    # ax.plot(z, c_inh(z), label = r'$c = \frac{1}{const. + f_i} -1 $' , c='darkgreen')

    ax.set_xlabel(r'$c$')
    # ax2.legend()
    # ticks = np.array([-np.pi, -np.pi/2, 0, np.pi/2 ,np.pi])
    # labels = [r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$']
    # ax.set_xticks(ticks)
    # ax.set_xticklabels(labels)

    cbar = plt.colorbar(cmap)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(r'sensitivity $s_\sigma$', rotation=270)

    fig.legend(bbox_to_anchor=(0.5, 0.9), loc='lower center')

    plt.show()


if __name__ == "__main__":
    main()
