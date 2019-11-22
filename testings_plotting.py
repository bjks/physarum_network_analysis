import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from matplotlib import rc
rc('text', usetex=True)

from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import curve_fit
from scipy import signal
from skimage.util.shape import view_as_windows

from scipy import LowLevelCallable

import numba
from numba import cfunc, carray
from numba.types import intc, CPointer, float64, intp, voidptr


import os
import itertools
import copy


from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *


from skimage import feature
from skimage import transform
from skimage.util.shape import view_as_windows
from skimage.feature import match_template


### plot ###
fig, axes = plt.subplots(2,1)
ax = axes.ravel()

shift = np.eye(10)
avg_shift = np.eye(10)*2

ax[0].set_title('phase shift')
ax[0].set_ylabel('space (pixel)')
ax[0].set_xlabel('time (s)')


ticks = np.array([-np.pi, -np.pi/2, 0, np.pi/2 ,np.pi])
labels = [r'$-\pi$',r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$']



im = ax[0].imshow(shift, origin='lower',
                aspect='auto', cmap='twilight_r')

im = ax[1].imshow(avg_shift, origin='lower',
                aspect='auto', cmap='twilight_r')

cbar = fig.colorbar(im, ax=ax.tolist(), orientation='vertical', ticks=ticks)
cbar.ax.set_yticklabels(labels)

plt.show()
