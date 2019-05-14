from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from skeleton import *
from correlation import *

import os
from scipy import signal
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import find_peaks


##########################################################################################

def fft_kymo(kymo):
    squared_kymo = crop_aligned_kymo(kymo)
    squared_kymo = detrend(squared_kymo, 40)
    kymo_f = np.fft.fft2(squared_kymo)
    # real = ndi.gaussian_filter(real, sigma=1)

    # real = np.fft.fftshift(kymo_f).real
    plt.imshow(np.transpose(squared_kymo))
    plt.show()

def power_spec(kymo1, kymo2, label, path_name):
    squared_kymo1_T = np.transpose(crop_aligned_kymo(kymo1))
    squared_kymo1_T = detrend(squared_kymo1_T, 40)

    squared_kymo2_T = np.transpose(crop_aligned_kymo(kymo2))
    squared_kymo2_T = detrend(squared_kymo2_T, 40)

    f1, Pxx_den1 = signal.periodogram(squared_kymo1_T)
    f2, Pxx_den2 = signal.periodogram(squared_kymo2_T)

    fig, axes = plt.subplots(2,2, figsize=(5,6))
    ax = axes.ravel()

    ax[0].imshow(Pxx_den1, vmin=0.03)
    ax[0].set_title('radius')

    ax[1].imshow(Pxx_den2, vmin=0.03)
    ax[1].set_title('concentration')


    ax[2].plot(f1, np.mean(Pxx_den1, axis=0))

    peaks, _ = find_peaks(np.mean(Pxx_den1, axis=0), distance=10, height=10)
    # ax[2].set_xticks(f1[peaks])
    ax[2].set_yticks([])
    ax[2].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[2].set_xlim([0,0.2])
    ax[2].set_xlabel('frequencies (1/frame)')
    ax[2].set_ylabel('power')

    ax[3].plot(f2, np.mean(Pxx_den2, axis=0))
    peaks, _ = find_peaks(np.mean(Pxx_den2, axis=0), distance=10)
    # ax[3].set_xticks(f2[peaks])
    ax[3].set_yticks([])
    ax[3].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[3].set_xlim([0,0.2])
    ax[3].set_xlabel('frequencies (1/frame)')

    plt.savefig(path_name + 'branch' + str(label) + '_power_spec.pdf')
    plt.close()




def main():
    set_keyword     = os.sys.argv[1]
    label           = os.sys.argv[2]
    color           = os.sys.argv[3]
    align_keyword   = 'reference_point'
    method          = 'inter_mean'
    max_offset      = 50

    set     = data(set_keyword, method=method, color=color)

    kymos = np.load(set.file_dat_set + '_branch_' + str(label) + '.npz')
    alignment=kymos['alignment']

    radii = align_kymo(kymos['kymograph_local_radii'], align_keyword, alignment=alignment)
    conce = align_kymo(kymos['kymograph_concentration'], align_keyword, alignment=alignment)

    path_name = set.file_plot_set + '_branch_' + str(label) + '/'

    fft_kymo(radii)
    power_spec(radii, conce, label, path_name)


if __name__ == '__main__':
    main()
