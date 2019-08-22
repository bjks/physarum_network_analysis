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
from analysis.plotting import *



############### Tools ###################
def mk_mising_dir(path_name):
    if not os.path.exists(path_name):
        os.mkdir(path_name)
    return path_name

############# FFT ############
def fft_kymo(signal, frame_int):

    kymo_f  = np.fft.rfft(signal)

    freq    = np.fft.rfftfreq(signal.shape[-1]) / frame_int
    fourier = np.mean(kymo_f.real, axis=0)
    plt.plot(freq, fourier)
    plt.show()


def bandpass(kymo, frame_int, min_freq=None, max_freq=None,
            min_period=None, max_period=None):

    # use period if given:
    if max_period != None:
        min_freq = 1./max_period

    if min_period != None:
        max_freq = 1./min_period

    kymo_f = np.fft.rfft(kymo)
    freq    = np.fft.rfftfreq(kymo.shape[-1], d=frame_int)

    if max_freq != None:
        kymo_f[:, (freq > max_freq)] = 0
    if min_freq != None:
        kymo_f[:, (freq < min_freq)] = 0

    return np.fft.irfft(kymo_f)


def dominant_freq(kymo, frame_int, min_freq=0.05, max_period=None):

    if max_period!=None:
        min_freq = 1./max_period

    f, Pxx_den_full = signal.periodogram(kymo, fs=1/frame_int)
    Pxx_den = np.mean(Pxx_den_full, axis=0)

    Pxx_den = Pxx_den[f>min_freq]
    f = f[f>min_freq]

    dominant_freq = f[np.argmax(Pxx_den)]
    return dominant_freq



def power_spec(kymo1, frame_int, file_name, title,
                min_freq=0.005, max_freq=0.05,          #periods between 200 and 20s
                min_period=None, max_period=None):

    if min_period!=None and max_period!=None:
        min_freq = 1./max_period
        max_freq = 1./min_period

    f1, Pxx_den1 = signal.periodogram(kymo1, fs=1/frame_int)
    Pxx_den1 = ndi.gaussian_filter(Pxx_den1, sigma=2)

    range = (f1 < max_freq) * (f1 > min_freq)

    Pxx_den1    = Pxx_den1[:,range]
    f1          = f1[range]


    fig, axes = plt.subplots(2,1, figsize=(5,6))
    ax = axes.ravel()

    ax[0].pcolormesh(Pxx_den1)
    ax[0].set_title('power spectrum ' + title)
    ax[0].set_xticks([])


    ax[1].plot(f1, np.mean(Pxx_den1, axis=0))
    # peaks, _ = find_peaks(np.mean(Pxx_den1, axis=0), distance=10)
    #
    # ax[1].set_xticks(f1[peaks])
    # ax[1].set_yticks([])
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax[1].set_xlabel('frequencies (1/s)')
    ax[1].set_ylabel('power')
    ax[1].set_xlim(f1[0], f1[-1])

    if SHOW:
        plt.show()

    plt.savefig(file_name + title + 'power.pdf', dpi=200)
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

    plt.savefig(file_name + title + 'hilbert.pdf', dpi=400)

    if SHOW:
        plt.show()
    plt.close()

    return instantaneous_phase, amplitude_envelope, instantaneous_frequency


###### Preparing Kymo ######
def get_kymo(kymos_data, keyword, frame_int, align_keyword='reference_point',
            min_period=40, max_period=130,
            times=(None,None), positions=(None, None)):

    alignment = kymos_data['alignment']

    kymo      = kymos_data[keyword]

    if align_keyword == None:
        return np.transpose(kymo)
                
    kymo = align_kymo(kymo, align_keyword, alignment=alignment)
    kymo = np.transpose(crop_aligned_kymo(kymo))

    kymo = kymo[positions[0]:positions[1], times[0]:times[1]]

    return kymo


def shift_radii(radii):
    return ndi.zoom(radii, (1,2), order=5)[1:]


def normalize_green(green, texas):
    green = ndi.zoom(green, (1,2), order=5)[:-1]
    texas = ndi.zoom(texas, (1,2), order=5)[1:]
    print(np.shape(texas))

    return green/texas




####################### Endo-ectoplasm problem ########################
def radius_dependence(radii, conce, frame_int):
    # g means 'global'
    radii_g = bandpass(radii, frame_int, min_period=500)
    show_im(radii_g)
    conce_g = bandpass(conce, frame_int, min_period=500)
    show_im(conce_g)

    # bin_edges = np.histogram_bin_edges(radii_g, bins='fd')
    bins = int(len(radii_g.flatten())/5000)

    c_of_R, edges_R,_   = stats.binned_statistic(radii_g.flatten(), conce_g.flatten(),
                                                    statistic='mean',
                                                    bins=bins)

    c_of_R = ndi.gaussian_filter(c_of_R, sigma=np.max(radii)/10)

    bin_R = edges_R[:-1] + np.diff(edges_R)/2

    plt.plot(bin_R, c_of_R)
    plt.xlabel("R")
    plt.ylabel("average concentration")
    plt.show()
    plt.close()

    # to account for extreme values of the unbandpassed(!) values:
    edges_R[0] = -np.inf
    edges_R[-1] = np.inf

    return edges_R, c_of_R


def substract_ecto_contribution(radii, concentration, frame_int):
    edges_R, c_of_R = radius_dependence(radii, concentration, frame_int)

    index = np.digitize(radii.ravel(), edges_R, right=True) - 1
    contribution = c_of_R[index].reshape(concentration.shape)
    show_im(contribution)

    return concentration - contribution


def plot_kymographs(kymos, titles, file_name, frame_int):

    for kymo, title in zip(kymos, titles):
        plt.title(title)
        plt.ylabel('space (pixel)')
        plt.xlabel('time (s)')
        plt.imshow(kymo)
        plt.colorbar()

        nx, nt = np.shape(kymo)[0], np.shape(kymo)[1]

        no_labels = 6 # how many labels to see on axis x
        step_t = int(nt / (no_labels - 1)) # step between consecutive labels
        t_positions = np.arange(0, nt, step_t) # pixel count at label position

        step_x = int(nx / (no_labels - 1)) # step between consecutive labels
        x_positions = np.arange(0, nx, step_x)

        plt.xticks(t_positions,  np.around(t_positions*frame_int, decimals=0))
        # plt.yticks(x_positions,  x_positions )

        plt.savefig(file_name + title + '_bandpasses_' + '.pdf', dpi=400)
        plt.close()


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



def sin_func(phase, A, phi):
    return A * np.cos(phase - phi)


def phase_average(phase, kymos, titles, colors, file_name):
    phase_shift = []
    for k, title, color in zip(kymos, titles, colors):
        bin_mean, bin_edges, bin_std = bin_acc_phase(phase, k, n_bins=15, norm=True)
        bin = bin_edges[:-1] + np.diff(bin_edges)/2

        phase_sample = np.linspace(bin[0], bin[-1], 100)
        popt, pcov   = curve_fit(sin_func, bin, bin_mean, p0=(-1,0))


        plt.plot(bin, bin_mean , label=title + ', shift: ' +
                    str( np.around(popt[-1], decimals=3) ), color=color)

        # plt.fill_between(bin, bin_mean - bin_std, bin_mean + bin_std, alpha=0.2, color=color)

        # min = bin[np.argmin(bin_mean)]
        # plt.axvline(x=min, linewidth=1, color=color)

        plt.plot(phase_sample, sin_func(phase_sample, *popt), '--', color=color)
        y1 = sin_func(phase_sample, popt[0], popt[1] + pcov[1,1]**0.5)
        y2 = sin_func(phase_sample, popt[0], popt[1] - pcov[1,1]**0.5)
        # plt.fill_between(phase_sample, y1, y2, color=color, alpha=0.15)

        plt.axvline(x=popt[-1], linewidth=1, color=color)
        phase_shift.append(popt[-1])

    plt.ylabel('average')
    plt.xlabel('phase')
    plt.legend()

    if SHOW:
        plt.show()
    plt.savefig(file_name + '_phase.pdf', dpi=400)
    plt.close()

    return phase_shift



############# CORRELATION #############
def gauss_detrend(kymo, r):
    global_structs = ndi.gaussian_filter(kymo, sigma=r)
    return kymo - global_structs


def correlate_phase(kymo1, kymo2, file_name, title, upsample_t=1, upsample_x=1, frame_int=1.0):

    kymo1 = ndi.zoom(kymo1, (upsample_x,upsample_t), order=5)
    kymo2 = ndi.zoom(kymo2, (upsample_x,upsample_t), order=5)
    # show_im(kymo1)

    image_product = np.fft.fft2(kymo1) * np.fft.fft2(kymo2).conj()
    cc_image = np.fft.ifft2(image_product)


    unshifted_correlation = cc_image.real
    correlation = np.fft.fftshift(cc_image).real

    nx, nt = np.shape(correlation)[0], np.shape(correlation)[1]
    max = np.argmax(unshifted_correlation[0])
    min = np.argmin(unshifted_correlation[0])
    print('Max: ', max/upsample_t, 'Min: ', min/upsample_t)
    if min > nt/2:
        min -= nt
    if max > nt/2:
        max -= nt


    #### plotting ####
    plt.imshow(correlation)

    plt.ylabel('space lag (pixel)')
    plt.xlabel('time lag (s)')


    max_in_s = np.around(max*frame_int/upsample_t, decimals=2)
    min_in_s = np.around(min*frame_int/upsample_t, decimals=2)

    plt.axvline(x=nt/2 + max, linewidth=0.5, color='k',
                label='Max: '+str(max_in_s) + ' s')
    plt.axvline(x=nt/2 + min, linewidth=0.5, color='r',
                label='Min: '+str(min_in_s) + ' s')

    no_labels = 6 # how many labels to see on axis x
    step_t = int(nt / (no_labels - 1)) # step between consecutive labels
    t_positions = np.arange(0, nt, step_t) # pixel count at label position

    step_x = int(nx / (no_labels - 1)) # step between consecutive labels
    x_positions = np.arange(0, nx, step_x)

    plt.xticks(t_positions,  (t_positions-nt/2)*frame_int/upsample_t)
    plt.yticks(x_positions,  (x_positions-nx/2)/upsample_x )

    plt.colorbar()
    plt.grid(linestyle='-', linewidth=0.4)
    plt.legend()

    if SHOW:
        plt.show()

    plt.tight_layout()
    plt.savefig(file_name + title + 'correlate.pdf', dpi=400)
    plt.close()


    print(min * frame_int/ upsample_t, max * frame_int/ upsample_t)
    return min * frame_int/ upsample_t, max * frame_int/ upsample_t




##########################################################################
################################ MAIN ####################################
##########################################################################
SHOW = False
SAVE = False

def main():
    set_keyword     = os.sys.argv[1].strip()
    color           = 'sep'

    align_keyword   = 'reference_point'
    method          = 'inter_mean'

    times           = None, None
    positions       = None, None

    substract_ecto  = False

    max_period      = 140
    range_freq      = 0.001, 0.003


    labels = range(len(data(set_keyword).seed_positions))
    for label in labels:
        set         = data(set_keyword, method=method, color=color)
        kymos_data  = np.load(set.file_dat_set + '_branch_' + str(label) + '.npz')

        path_name = mk_mising_dir(set.file_plot_set + '_branch_' + str(label) + '/')
        file_name = path_name + '/branch_' + str(label) + '_'

        ##################################################################
        ########################## Collect data ##########################
        ##################################################################
        radii = get_kymo(kymos_data, 'kymo_local_radii', set.frame_int,
                        align_keyword, times=times, positions=positions)

        kymo_c_green    = get_kymo(kymos_data, 'kymo_c_green', set.frame_int,
                        align_keyword, times=times, positions=positions)
        kymo_inner_green = get_kymo(kymos_data, 'kymo_inner_green',  set.frame_int,
                        align_keyword, times=times, positions=positions)
        kymo_outer_green = get_kymo(kymos_data, 'kymo_outer_green',  set.frame_int,
                        align_keyword, times=times, positions=positions)

        kymo_c_texas    = get_kymo(kymos_data, 'kymo_c_texas', set.frame_int,
                        align_keyword, times=times, positions=positions)
        kymo_inner_texas = get_kymo(kymos_data, 'kymo_inner_texas',  set.frame_int,
                        align_keyword, times=times, positions=positions)
        kymo_outer_texas = get_kymo(kymos_data, 'kymo_outer_texas',  set.frame_int,
                        align_keyword, times=times, positions=positions)


        titles = ['radius', 'green', 'texas']
        kymos  = [radii, kymo_c_green, kymo_c_texas]

        plot_kymographs(kymos, titles, file_name, set.frame_int)

        radii = shift_radii(radii)
        conce = normalize_green(kymo_c_green, kymo_c_texas)
        inner = normalize_green(kymo_inner_green, kymo_inner_texas)
        outer = normalize_green(kymo_outer_green, kymo_outer_texas)
        set.frame_int/=2.

        print(set.frame_int)

        freq = dominant_freq(radii, set.frame_int, max_period=max_period)
        print('Dominant frequency: ', freq)


        ####################################################
        ########### substract  Ca_global(R) ################
        ####################################################
        if substract_ecto:
            show_im(conce)
            conce = substract_ecto_contribution(radii, conce, set.frame_int)
            inner = substract_ecto_contribution(radii, inner, set.frame_int)
            outer = substract_ecto_contribution(radii, outer, set.frame_int)
            show_im(conce)
        ##################################################################
        ################# calc phase of oscillation ######################
        ##################################################################
        radius_base = bandpass(radii, set.frame_int,
                                min_freq=freq-range_freq[0],
                                max_freq=freq+range_freq[1])

        conce_base = bandpass(conce, set.frame_int,
                                min_freq=freq-range_freq[0],
                                max_freq=freq+range_freq[1])

        phase_radius, _, _ = extract_phase(radius_base, file_name, 'radius')
        phase_conce, _, _  = extract_phase(conce_base, file_name, 'concentration')



        ##################################################################
        ############ bandpass kymographs before binning ##################
        ##################################################################
        radii = bandpass(radii, set.frame_int,  min_freq=freq-range_freq[0],
                                                max_freq=freq+range_freq[1])
        conce = bandpass(conce, set.frame_int,  min_freq=freq-range_freq[0],
                                                max_freq=freq+range_freq[1])
        inner = bandpass(inner, set.frame_int,  min_freq=freq-range_freq[0],
                                                max_freq=freq+range_freq[1])
        outer = bandpass(outer, set.frame_int,  min_freq=freq-range_freq[0],
                                                max_freq=freq+range_freq[1])


        ##################################################################
        ######################## Phase dependence ########################
        ##################################################################
        titles = ['radius', 'concentration', 'inner', 'outer']
        colors = ['orange', 'blue', 'darkblue', 'lightskyblue']
        kymos  = [radii, conce, inner, outer]

        plot_kymographs(kymos, titles, file_name, set.frame_int)


        phase_shift = phase_average(phase_radius, kymos, titles, colors,
                                    file_name)


        ##################################################################
        #################### Fourier spectrum  ###########################
        ##################################################################
        power_spec(radii, set.frame_int, file_name, 'radius')
        power_spec(conce, set.frame_int, file_name, 'concentration')

        ##################################################################
        ######################### correlations ###########################
        ##################################################################
        min_k, max_k = correlate_phase(radii, conce,
                                        file_name, 'kymos',
                                        upsample_t=10, upsample_x=3,
                                        frame_int=set.frame_int)


        min_p, max_p = correlate_phase(phase_radius, phase_conce,
                                        file_name, 'phases',
                                        upsample_t=10, upsample_x=3,
                                        frame_int=set.frame_int)


        ##################################################################
        ######################### Save in txt ############################
        ##################################################################
        if SAVE:
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
