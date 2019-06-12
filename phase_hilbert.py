import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import curve_fit

import os

from scipy import signal
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import find_peaks

from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *


############# FFT ############
def fft_kymo(signal, frame_int):

    kymo_f  = np.fft.rfft(signal)

    freq    = np.fft.rfftfreq(signal.shape[-1]) / frame_int
    fourier = np.mean(kymo_f.real, axis=0)
    plt.plot(freq, fourier)
    plt.show()


def bandpass(signal, frame_int,
            lower_cutoff=0.005, upper_cutoff=0.08): #periods between 20 and 200s

    kymo_f = np.fft.rfft(signal)
    freq    = np.fft.rfftfreq(signal.shape[-1], d=frame_int)

    kymo_f[:, (freq>upper_cutoff)] = 0
    kymo_f[:, (freq<lower_cutoff)] = 0

    return np.fft.irfft(kymo_f)


def power_spec(kymo1, frame_int, file_name, title, max_freq=0.05):

    f1, Pxx_den1 = signal.periodogram(kymo1, fs=1/frame_int)
    range = f1 < max_freq

    Pxx_den1    = Pxx_den1[:,range]
    f1          = f1[range]


    fig, axes = plt.subplots(2,1, figsize=(5,6))
    ax = axes.ravel()

    ax[0].pcolormesh(Pxx_den1)
    ax[0].set_title('power spectrum ' + title)
    ax[0].set_xticks([])


    ax[1].plot(f1, np.mean(Pxx_den1, axis=0))
    peaks, _ = find_peaks(np.mean(Pxx_den1, axis=0), distance=10)

    ax[1].set_xticks(f1[peaks])
    # ax[1].set_yticks([])
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    ax[1].set_xlabel('frequencies (1/s)')
    ax[1].set_ylabel('power')
    ax[1].set_xlim(f1[0], f1[-1])

    if SHOW:
        plt.show()

    plt.savefig(file_name + title + 'power.pdf')
    plt.close()

############## HILBERT ###############
def extract_phase(signal, file_name, title):

    analytic_signal = hilbert(signal)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.angle(analytic_signal)
    # instantaneous_phase = np.unwrap(instantaneous_phase, axis=0)
    instantaneous_frequency = np.diff(instantaneous_phase)

    fig, axes = plt.subplots(1,4, figsize=(6, 3), sharex=True, sharey=True)
    ax = axes.ravel()
    ax[0].set_title(title)
    c = ax[0].imshow(signal)

    ax[1].set_title('amplitude')
    ax[1].imshow(amplitude_envelope)

    ax[2].set_title('phase')
    ax[2].imshow(instantaneous_phase)

    ax[3].set_title('frequence')
    ax[3].imshow(instantaneous_frequency)

    plt.savefig(file_name + title + 'hilbert.pdf')

    if SHOW:
        plt.show()
    plt.close()

    return instantaneous_phase, amplitude_envelope, instantaneous_frequency


###### Preparing Kymo Array ######
def kymo_prep(kymo, align_keyword=None, alignment=None):

    if align_keyword == None:
        return np.transpose(kymo)

    aligned_kymo = align_kymo(kymo, align_keyword, alignment=alignment)
    return np.transpose(crop_aligned_kymo(aligned_kymo))


def collect_kymos(kymos_data, keywords, frame_int, align_keyword='reference_point'):
    kymos = []
    alignment       = kymos_data['alignment']

    for k in keywords:
        kymo = kymos_data[k]
        kymo = kymo_prep(kymo, align_keyword, alignment)
        kymo = bandpass(kymo, frame_int)
        # fft_kymo(kymo, set.frame_int)
        # show_im(kymo)
        kymos.append(kymo)

    return kymos


def gauss_detrend(kymo, r):
    global_structs = ndi.gaussian_filter(kymo, sigma=r)
    return kymo - global_structs


#################### BINNING ####################
def bin_acc_phase(phase, x, n_bins=10, norm=True):
    phase_f = phase.flatten()
    x_f     = x.flatten()

    x_mean, bin_edges,_   = stats.binned_statistic(phase_f, x_f,
                                                    statistic='mean',
                                                    bins=n_bins)
    x_std,_, _            = stats.binned_statistic(phase_f, x_f,
                                                    statistic=np.std,
                                                    bins=n_bins)

    if norm:
        x_max = np.max(x_mean)
        x_mean /= x_max
        x_std /= x_max

    return x_mean, bin_edges, x_std


def phase_average(phase, kymos, titles, colors, file_name):
    phase_shift = []
    for k, title, color in zip(kymos, titles, colors):
        bin_mean, bin_edges, bin_std = bin_acc_phase(phase, k, n_bins=20)
        bin = bin_edges[:-1] + np.diff(bin_edges)/2

        phase_sample = np.linspace(bin[0], bin[-1], 100)
        popt, _ = curve_fit(sin_func, bin, bin_mean)

        plt.plot(phase_sample, sin_func(phase_sample, *popt), '--', color=color)

        plt.plot(bin, bin_mean , label=title + ', shift: ' +
                    str( np.around(popt[-1], decimals=3) ),
                    color=color)
        # plt.fill_between(bin, r_mean - r_std, r_mean + r_std, alpha=0.2, color='orange')

        plt.axvline(x=popt[-1], linewidth=1, color=color)
        phase_shift.append(popt[-1])

    plt.ylabel('average (a.u.)')
    plt.xlabel('phase')
    plt.legend()

    if SHOW:
        plt.show()
    plt.savefig(file_name + '_phase.pdf')
    plt.close()


    return phase_shift



############# CORRELATION #############
def correlate_phase(kymo1, kymo2, file_name, title, upsample=1, frame_int=1.0):

    kymo1 = ndi.zoom(kymo1, upsample, order=5)
    kymo2 = ndi.zoom(kymo2, upsample, order=5)

    image_product = np.fft.fft2(kymo1) * np.fft.fft2(kymo2).conj()
    cc_image = np.fft.ifft2(image_product)


    unshifted_correlation = cc_image.real
    correlation = np.fft.fftshift(cc_image).real

    nx, nt = np.shape(correlation)[0], np.shape(correlation)[1]
    max = np.argmax(unshifted_correlation[0])
    min = np.argmin(unshifted_correlation[0])
    print('Max: ', max/upsample, 'Min: ', min/upsample)
    if min > nt/2:
        min -= nt
    if max > nt/2:
        max -= nt


    #### plotting ####
    plt.imshow(correlation)

    plt.ylabel('space lag (frame)')
    plt.xlabel('time lag (s)')


    plt.axvline(x=nt/2 + max, linewidth=0.5, color='k',
                label='Max: '+str(max*frame_int/upsample) + ' s')
    plt.axvline(x=nt/2 + min, linewidth=0.5, color='r',
                label='Min: '+str(min*frame_int/upsample) + ' s')

    no_labels = 5 # how many labels to see on axis x
    step_t = int(nt / (no_labels - 1)) # step between consecutive labels
    t_positions = np.arange(0, nt, step_t) # pixel count at label position

    step_x = int(nx / (no_labels - 1)) # step between consecutive labels
    x_positions = np.arange(0, nx, step_x)

    plt.xticks(t_positions,  (t_positions-nt/2)*frame_int/upsample)
    plt.yticks(x_positions,  (x_positions-nx/2)/upsample )

    plt.colorbar()
    plt.grid(linestyle='-', linewidth=0.4)
    plt.legend()

    if SHOW:
        plt.show()

    plt.savefig(file_name + title + 'correlate.pdf')
    plt.close()

    return min * frame_int/ upsample, max * frame_int/ upsample



def sin_func(phase, A, phi):
    return A * np.cos(phase - phi)



##########################################################################
################################ MAIN ####################################
##########################################################################
SHOW = False

def main():
    set_keyword     = os.sys.argv[1]
    colors           = os.sys.argv[2:]

    align_keyword   = 'reference_point'
    method          = 'inter_mean'



    labels = range(len(data(set_keyword).seed_positions))
    for color in colors:
        for label in labels:
            set         = data(set_keyword, method=method, color=color)
            kymos_data  = np.load(set.file_dat_set + '_branch_' + str(label) + '.npz')


            keywords = ['kymograph_local_radii',
                        'kymograph_concentration',
                        'kymograph_inner',
                        'kymograph_outer']

            titles = ['radius', 'concentration', 'inner', 'outer']
            colors = ['orange', 'blue', 'darkblue', 'lightskyblue']

            kymos = collect_kymos(kymos_data, keywords, set.frame_int, align_keyword)

            path_name = set.file_plot_set + '_branch_' + str(label) + '/'
            if not os.path.exists(path_name):
                os.mkdir(path_name)
            file_name = path_name + '/branch_' + str(label) + '_'

            ########################## analysis #############################
            ###### power spec   ######
            power_spec(kymos[0], set.frame_int, file_name, 'radius')
            power_spec(kymos[1], set.frame_int, file_name, 'concentration')

            ###### phases       ######
            phase, amp, freq        = extract_phase(kymos[0], file_name, 'radius')
            phase_c, amp_c, freq_c  = extract_phase(kymos[1], file_name, 'concentration')

            phase_shift = phase_average(phase, kymos, titles, colors, file_name)

            ###### correlations ######
            min_k, max_k = correlate_phase(kymos[0], kymos[1], file_name, 'kymos',
                                            upsample=5, frame_int=set.frame_int)
            min_p, max_p = correlate_phase(phase, phase_c, file_name,
                                            'phases', upsample=2)


            ######################### Save in txt ###########################
            if not os.path.exists('time_shifts_data_sets'):
                os.mkdir('time_shifts_data_sets')

            data_sets_summary = 'time_shifts_data_sets/time_shift_' + color +'.txt'
            if not os.path.isfile(data_sets_summary):
                with open(data_sets_summary, "w") as out_var:
                    out_var.write('# data_set \t label \t phase_radius' +
                                    '\t phase_conc \t phase_inner' +
                                    '\t phase_outer \t min_kymo \t max_kymo' +
                                     '\t min_phase \t max_phase \n')

            with open(data_sets_summary, "a") as out_var:
                out_var.write(set_keyword +
                                '\t' + str(label) +
                                '\t' + str(phase_shift[0]) +
                                '\t' + str(phase_shift[1]) +
                                '\t' + str(phase_shift[2]) +
                                '\t' + str(phase_shift[3]) +
                                '\t' + str(min_k) +
                                '\t' + str(max_k) +
                                '\t' + str(min_p) +
                                '\t' + str(max_p) + '\n')


if __name__ == '__main__':
    main()
