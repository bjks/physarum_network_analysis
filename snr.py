import matplotlib.pyplot as plt
from analysis.network_analysis import *
from analysis.data_sets import *
from scipy import stats
import os


def snr(data_sets, color):

    masks    = [np.load(set.file_dat + '.npz')['mask'].astype(bool) for set in data_sets]

    if color=='texas':
        raw = [read_file(set.file_raw2) for set in data_sets]
    elif color=='green':
        raw = [read_file(set.file_raw1) for set in data_sets]

    snr=[]

    for r, m in zip(raw, masks):
        signal  = np.true_divide(np.sum(r*m), np.sum(m))
        back    = np.true_divide(np.sum(r*np.invert(m)), np.sum(np.invert(m)))
        snr.append(signal/back)
        print(2)

    return snr





def main():
    set_keyword = os.sys.argv[1]
    step        = int(os.sys.argv[2])
    method      = 'inter_mean'

    # set_inds = np.arange(data(set_keyword).first, data(set_keyword).last, step)
    set_inds = np.arange(data(set_keyword).first, 10, step)


    data_sets = [data(set_keyword, i, method=method) for i in set_inds]


    dt = data(set_keyword).frame_int
    time = [i * dt for i in range(len(set_inds))]

    snr_g = snr(data_sets, 'green')
    snr_t = snr(data_sets, 'texas')


    plt.plot(time, snr_g, color = 'green', label = 'SNR green')
    plt.plot(time, snr_t, color = 'red', label = 'SNR texas red')

    plt.ylabel('SNR')
    plt.xlabel('time (s)')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
