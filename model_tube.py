from skimage.feature import canny

import numpy as np
import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
from correlation import *
from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter

from multiprocessing.dummy import Pool as ThreadPool

def wall_height(r, d):
    R1 = 1 - d
    hw = np.where(r < R1, np.sqrt(1 - r**2) - np.sqrt(R1**2 - r**2), np.sqrt(1 - r**2))
    return hw

def tube_height(R, r):
    return np.sqrt(R**2 - r**2)

def main():

    r = np.linspace(0, 1, 1000)
    r = r[:-1]
    plt.plot(r, tube_height(1, r), c='red')
    plt.plot(r, wall_height(r, 0.1), c='green')
    print(np.mean(tube_height(2, 2*r))/ np.mean(tube_height(1, r)))
    plt.plot(r, wall_height(r, 0.1)/tube_height(1, r), c='blue')

    plt.savefig('model_tube.pdf')

    plt.show()

if __name__ == '__main__':
    main()
