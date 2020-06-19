import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import itertools
import os

import sys
sys.path.append("..")

from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.kymograph_analysis import *
from analysis.plotting import *

def plot_bf(file, scaling, cmap, thick=False, cbar=True):

    image = read_file(file)[:800,:]

    fig, ax = plt.subplots(frameon=False)

    im = ax.imshow(image, cmap=cmap)

    # correction of scalebar_size to get 1,10,50... um
    # used in pretty mode
    allowed_scalebars = np.append([10**i for i in range(10)],
                                [5* 10**i for i in range(10)] )
    scalebar_size = closest_number(allowed_scalebars,
                                    400 * scaling)
    pixels_scalebar = scalebar_size/scaling
    scalebar = AnchoredSizeBar(ax.transData,
                       pixels_scalebar, str(scalebar_size) + r' $\mu$m',
                       'lower right',
                       pad=0.3,
                       color='black',
                       frameon=False,
                       size_vertical=1)


    ax.add_artist(scalebar)
    ax.axis('off')

    ax.set_yticks([])
    ax.set_xticks([])
    # ax.set_title(title)

    # cbar = fig.colorbar(im)
    fig.tight_layout()

    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.set_ylabel('intensity (a.u.)', rotation=270)

    savefile = file + 'scalebar' + '.pdf'
    print(savefile)
    plt.savefig(savefile, dpi=600, bbox_inches='tight')
    plt.close()



def main():
    file = os.sys.argv[1].strip()
    scaling = float(os.sys.argv[2].strip())
    plot_bf(file, scaling, 'Greys_r', cbar=False)

if __name__ == '__main__':
    main()
