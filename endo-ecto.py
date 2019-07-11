import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.data_sets import *
from scipy import stats
import os


def diff(texas1, texas2):
    texas1 = ndi.gaussian_filter(texas1, sigma=2)
    texas2 = ndi.gaussian_filter(texas2, sigma=2)

    # diff = np.true_divide(np.absolute(texas2-texas1), (texas1+texas2))

    diff= np.absolute(texas2-texas1)
    return diff * texas1.astype(bool) * texas2.astype(bool)

def main():
    set_keyword = os.sys.argv[1]
    step        = int(os.sys.argv[2])
    method      = 'inter_mean'

    # set_inds = np.arange(data(set_keyword).first, data(set_keyword).last, step)
    set_inds = np.arange(data(set_keyword).first, 10, step)


    data_sets = [data(set_keyword, i, method=method) for i in set_inds]

    mask    = [np.load(set.file_dat + '.npz')['mask'] for set in data_sets]
    texas   = [np.load(set.file_dat + '.npz')['texas_clean'] for set in data_sets]

    for i in range(len(set_inds-1)):
        plt.imshow(diff(texas[i], texas[i+1]))
        plt.show()


if __name__ == '__main__':
    main()
