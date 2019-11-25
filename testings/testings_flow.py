import numpy as np
import matplotlib.pyplot as plt


import sys
sys.path.append("..")
# from analysis.network_analysis import *
# from analysis.skeleton_analysis import *
# from analysis.data_sets import *
from analysis.flow_analysis import *

from timeit import default_timer as timer
# import os
# from scipy.ndimage.filters import generic_filter
# from phase_hilbert import *

# from multiprocessing.dummy import Pool as ThreadPool
###############



def main():

    im = np.clip(np.random.randn(250, 250), -5, 1)

    frame_a = im[6:]
    frame_b = im[3:-3]
    frame_c = im[:-6]



    mask = np.ones_like(frame_a)
    mask[:,-100:] = 0


    start_t = timer()
    x, y = average_flow_over_frames([frame_a, frame_b, frame_c], mask, upsample=1,
                            window_size=20, sampling=20, search_extend=20)

    print('time: ', timer()-start_t, ' s')



    plt.imshow(frame_a, origin='lower')
    plt.show()

    plt.imshow(frame_b, origin='lower')
    plt.show()

    plt.imshow(frame_c, origin='lower')
    plt.show()

    plt.imshow(mask, origin='lower')
    plt.show()

    show_flow(x,y)



if __name__ == '__main__':
    main()
