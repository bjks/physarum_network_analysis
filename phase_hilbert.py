from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *

from scipy.signal import hilbert


from skimage.feature import register_translation
import os

def phase2d(intensity):
    h = scipy.signal.hilbert2(intensity)
    # analytical = intensity + np.j * h
    return h

def main():
    set_keyword     = os.sys.argv[1]
    color           = os.sys.argv[2]

    align_keyword   = 'reference_point'
    method          = 'inter_mean'
    data_sets = [data(set_keyword, i, method, color=color) for i in range(data(set_keyword).first, data(set_keyword).last)]

    local_radii   = np.load(set.file_dat + '.npz')['local_radii']

    h = phase2d(local_radii)
    plt.imshow(h)
    plt.show()




if __name__ == '__main__':
    main()
