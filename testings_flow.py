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
    # im = invert_bf(read_file('../image_data/2019-09-13/TR_CG1_flow_pgol2/TR_CG1_flow_pgol2_t001c2.tif'))
    # mask = create_mask(im, 10)

    # show_im(mask)
    # show_im(im)

    im = np.eye(100,100)

    frame_a = im[6:]
    frame_b = im[3:-3]
    frame_c = im[:-6]


    # show_im(frame_a)
    # show_im(frame_b)
    # show_im(frame_c)


    mask = np.ones_like(frame_a)

    x, y = average_flow_over_frames([frame_a, frame_b, frame_c], mask, upsample=1,
                            window_size=20, sampling=50, search_extend=20,
                            corr_mode='valid')

    plt.imshow(frame_a)
    plt.show()
    dy = np.arange(np.shape(frame_a)[0])
    dx = np.arange(np.shape(frame_a)[1])
    plt.quiver(x,y, scale=0.1)

    plt.show()


if __name__ == '__main__':
    main()
