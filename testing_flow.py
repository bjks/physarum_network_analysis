from skimage.feature import canny

import numpy as np
import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.data_sets import *
from analysis.flow_analysis import *

from timeit import default_timer as timer
import os
from scipy.ndimage.filters import generic_filter
from phase_hilbert import *

from multiprocessing.dummy import Pool as ThreadPool
###############



def main():
    im = (np.random.rand(2000, 2000) * 10).astype(int)

    frame_a = im[2:]
    frame_b = im[1:-1]
    frame_c = im[:-2]


    mask = np.ones_like(frame_b)
    mask[:500] = 0
    print(mask)

    x, y = average_flow_over_frames([frame_a, frame_b, frame_c], mask, upsample=1,
                            window_size=40, sampling=20, search_extend=5)

    plt.quiver(x,y)
    plt.show()

if __name__ == '__main__':
    main()
