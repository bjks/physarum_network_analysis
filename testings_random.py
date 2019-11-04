
import numpy as np
import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter
from phase_hilbert import *

from multiprocessing.dummy import Pool as ThreadPool

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from numpy.random import randn

from skimage.feature import match_template


def main():

    fig, ax = plt.subplots()

    data = np.clip(randn(250, 250), -5, 1)


    corr = match_template(data, data[10:-10, 20:-20])
    cax = ax.imshow(corr, interpolation='nearest', cmap=cm.coolwarm, aspect='auto')


    ax.axvline(x=0.2, linewidth=0.5, color='r',)

    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    ticks = np.append( np.arange(-5, 5, 0.5), np.mean(data)  )
    labels = ticks.astype(str)
    labels[-1]='    mean'
    cbar = fig.colorbar(cax, ticks=ticks)
    cbar.ax.set_yticklabels(labels)  # vertically oriented colorbar


    plt.show()

    plt.plot(corr[10])
    plt.show()


if __name__ == '__main__':
    main()
