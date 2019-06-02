import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, chirp
from scipy import ndimage as ndi
import os

from scipy import signal
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import find_peaks

from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *



def extract_phase(signal, phase=False):

    analytic_signal = hilbert(signal)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.angle(analytic_signal)
    # instantaneous_phase = np.unwrap(instantaneous_phase, axis=0)
    instantaneous_frequency = np.diff(instantaneous_phase)
    if not phase:
        return instantaneous_phase, amplitude_envelope, instantaneous_frequency

    return instantaneous_phase


def fft_kymo(signal):
    kymo_f  = np.fft.rfft(signal)

    freq    = np.fft.rfftfreq(signal.shape[-1])
    fourier = np.mean(kymo_f.real, axis=0)
    plt.plot(freq, fourier)
    plt.show()

def bandpass(signal):
    kymo_f = np.fft.rfft(signal)

    kymo_f[:,   :10]     =0
    kymo_f[:,   -10:]   =0

    return np.fft.irfft(kymo_f)


def power_spec(kymo1, kymo2):

    f1, Pxx_den1 = signal.periodogram(kymo1)
    f2, Pxx_den2 = signal.periodogram(kymo2)

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
    plt.show()

    # plt.savefig(path_name + 'branch' + str(label) + '_power_spec.pdf')
    plt.close()

def kymo_prep(kymo, align_keyword=None, alignment=None):
    if align_keyword == None:
        return np.transpose(kymo)

    aligned_kymo = align_kymo(kymo, align_keyword, alignment=alignment)
    return np.transpose(crop_aligned_kymo(aligned_kymo))

def correlate_phase(kymo1, kymo2, upsample=1):
    kymo1 = ndi.zoom(kymo1, upsample, order=5)
    kymo2 = ndi.zoom(kymo2, upsample, order=5)

    image_product = np.fft.fft2(kymo1) * np.fft.fft2(kymo2).conj()
    cc_image = np.fft.ifft2(image_product)


    unshifted_correlation = cc_image.real
    correlation = np.fft.fftshift(cc_image).real

    x_range, t_range = np.shape(correlation)[0], np.shape(correlation)[1]
    max = np.argmax(unshifted_correlation[0])
    min = np.argmin(unshifted_correlation[0])
    print('Max: ', max/upsample, 'Min: ', min/upsample)
    if min > t_range/2:
        min -= t_range
    if max > t_range/2:
        max -= t_range

    plt.imshow(correlation)

    plt.ylabel('space lag')
    plt.xlabel('time lag')

    plt.axvline(x=t_range/2 + max, linewidth=0.5, color='k', label='Max: '+str(max/upsample))
    plt.axvline(x=t_range/2 + min, linewidth=0.5, color='r', label='Min: '+str(min/upsample))

    plt.xticks([t_range/2], (0,) )
    plt.yticks([x_range/2], (0,) )

    plt.colorbar()
    plt.grid(linestyle='-', linewidth=0.4)
    plt.legend()
    plt.show()

def main():

    set_keyword     = os.sys.argv[1]
    colors           = os.sys.argv[2:]

    align_keyword   = 'reference_point'
    method          = 'inter_mean'



    labels = range(len(data(set_keyword).seed_positions))
    for order in colors:
        for label in labels:
            set     = data(set_keyword, method=method, color=order)
            kymos = np.load(set.file_dat_set + '_branch_' + str(label) + '.npz')

            radii           = kymos['kymograph_local_radii']
            concentration   = kymos['kymograph_concentration']
            alignment       = kymos['alignment']

            radii           = kymo_prep(radii, align_keyword, alignment)
            radii = bandpass(radii)
            fft_kymo(radii)

            concentration   = kymo_prep(concentration, align_keyword, alignment)
            concentration = bandpass(concentration)
            fft_kymo(concentration)


            phase1, amp1, freq1 = extract_phase(radii, phase=False)
            phase2, amp2, freq2 = extract_phase(concentration, phase=False)

            power_spec(radii, concentration)

#########
            fig, axes = plt.subplots(2,4, figsize=(6, 6), sharex=True, sharey=True)
            ax = axes.ravel()
            ax[0].set_title('radii')
            c = ax[0].imshow(radii)

            ax[1].set_title('amplitude')
            ax[1].imshow(amp1)

            ax[2].set_title('phase freq')
            ax[2].imshow(phase1)

            ax[3].set_title('frequence')
            ax[3].imshow(freq1)
#########
            ax[4].set_title('concentration')
            ax[4].imshow(concentration)

            ax[5].imshow(amp2)

            ax[6].imshow(phase2)

            ax[7].imshow(freq2)

            plt.show()

            correlate_phase(phase1, phase2)


if __name__ == '__main__':
    main()
