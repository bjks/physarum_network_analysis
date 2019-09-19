import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import curve_fit

import os

from scipy import signal
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
                min_period=None, max_period=None,
                mark_freq=None, logscale=False):

    if min_period!=None and max_period!=None:
        min_freq = 1./max_period
        max_freq = 1./min_period

    f1, Pxx_den1 = signal.periodogram(kymo1, fs=1./frame_int)
    # Pxx_den1 = ndi.gaussian_filter(Pxx_den1, sigma=2)


    range = (f1 < max_freq) * (f1 > min_freq)

    Pxx_den1    = Pxx_den1[:,range]
    f1          = f1[range] * 1e3 ### go to mHz

    if logscale:
        Pxx_den1 = np.log(Pxx_den1)

    fig, axes = plt.subplots(2,1, figsize=(5,6))
    ax = axes.ravel()

    im = ax[0].imshow(Pxx_den1, aspect='auto')
    ax[0].set_title('power spectrum ' + title)
    ax[0].set_xticks([])


    divider = make_axes_locatable(ax[0])
    fig.colorbar(im, cax = divider.append_axes('right',size='5%', pad=0.05),
                 orientation='vertical')



    ax[1].plot(f1, np.mean(Pxx_den1, axis=0))
    # peaks, _ = find_peaks(np.mean(Pxx_den1, axis=0), distance=10)
    #
    # ax[1].set_xticks(f1[peaks])
    # ax[1].set_yticks([])
    if mark_freq != None:
        mark_freq *= 1e3
        ax[1].axvline(x=mark_freq, linewidth=1, color='red',
            label='dominant frequency: ' +
            str(np.around(mark_freq, decimals=2)) + ' mHz')

        ax[1].legend()

    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].set_xlabel('frequency (mHz)')
    if logscale:
        ax[1].set_ylabel('log(power)')
    else:
        ax[1].set_ylabel('power')

    ax[1].set_xlim(f1[0], f1[-1])
    fig.tight_layout()

    if SHOW:
        plt.show()

    plt.savefig(file_name + title + 'power.pdf', dpi=200)
    plt.close()


############## HILBERT ###############
def extract_phase(signal, file_name, title, frame_int):

    analytic_signal = hilbert(signal)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.angle(analytic_signal)
    # instantaneous_phase = np.unwrap(instantaneous_phase, axis=1)
    instantaneous_frequency = np.diff(instantaneous_phase)/(2.0*np.pi)*frame_int

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
def get_kymo(kymos_data, keyword, align_keyword='reference_point',
            times=(None,None), positions=(None, None), file_name=None, title=None, plot=False):

    alignment = kymos_data['alignment']

    kymo      = kymos_data[keyword]

    if align_keyword == None:
        return np.transpose(kymo)

    kymo = align_kymo(kymo, align_keyword, alignment=alignment)

    if plot:
        fig, ax = plt.subplots()

        ax.set_title(title)
        ax.set_ylabel('space (pixel)')
        ax.set_xlabel('time (frame)')
        im = ax.imshow(np.transpose(kymo))
        divider = make_axes_locatable(ax)
        fig.colorbar(im, cax = divider.append_axes('right', size='5%', pad=0.05),
                     orientation='vertical')
        fig.tight_layout()

        plt.savefig(file_name + title + '.pdf', dpi=400)
        plt.close()

    kymo = np.transpose(crop_aligned_kymo(kymo))
    kymo = kymo[positions[0]:positions[1], times[0]:times[1]]

    return kymo


def shift_radii(radii, symm_setup, frame_int):
    if symm_setup:
        return ndi.zoom(radii, (1,2), order=5)[1:], frame_int/2.
    else:
        return radii, frame_int

def normalize_green(green, texas, symm_setup):
    if symm_setup:
        green = ndi.zoom(green, (1,2), order=5)[:-1]
        texas = ndi.zoom(texas, (1,2), order=5)[1:]

    else:
        pass

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

        # nx, nt = np.shape(kymo)[0], np.shape(kymo)[1]
        #
        # no_labels = 6 # how many labels to see on axis x
        # step_t = int(nt / (no_labels - 1)) # step between consecutive labels
        # t_positions = np.arange(0, nt, step_t) # pixel count at label position
        #
        # step_x = int(nx / (no_labels - 1)) # step between consecutive labels
        # x_positions = np.arange(0, nx, step_x)
        # ax.set_xticks(t_positions, np.around(t_positions*frame_int, decimals=0))
        # plt.yticks(x_positions,  x_positions )


        fig, ax = plt.subplots()
        im = ax.imshow(kymo)

        ax.set_xlim(left=0, right=np.shape(kymo)[1])
        locs = ax.get_xticks()
        ax.set_xticks(locs)
        ax.set_xticklabels(np.around(locs * frame_int, decimals=0).astype(int))

        ax.set_title(title)
        ax.set_ylabel('space (pixel)')
        ax.set_xlabel('time (s)')

        divider = make_axes_locatable(ax)
        fig.colorbar(im, cax = divider.append_axes('right', size='5%', pad=0.05),
                     orientation='vertical')

        plt.savefig(file_name + title + '.pdf', dpi=400)
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

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    phase_shift = []

    for k, title, color in zip(kymos, titles, colors):
        bin_mean, bin_edges, bin_std = bin_acc_phase(phase, k, n_bins=15, norm=False)
        bin = bin_edges[:-1] + np.diff(bin_edges)/2


        phase_sample = np.linspace(bin[0], bin[-1], 100)
        popt, pcov   = curve_fit(sin_func, bin, bin_mean, p0=(-1,0))

        if title == titles[0]:
            axis = ax1
        else:
            axis = ax2

        axis.plot(bin, bin_mean , label=title + ', shift: ' +
                    str( np.around(popt[-1], decimals=3) ), color=color)

        # plt.fill_between(bin, bin_mean - bin_std, bin_mean + bin_std, alpha=0.2, color=color)

        # min = bin[np.argmin(bin_mean)]
        # plt.axvline(x=min, linewidth=1, color=color)

        axis.plot(phase_sample, sin_func(phase_sample, *popt), '--', color=color)
        y1 = sin_func(phase_sample, popt[0], popt[1] + pcov[1,1]**0.5)
        y2 = sin_func(phase_sample, popt[0], popt[1] - pcov[1,1]**0.5)
        # plt.fill_between(phase_sample, y1, y2, color=color, alpha=0.15)

        axis.axvline(x=popt[-1], linewidth=1, color=color)
        phase_shift.append(popt[-1])


    ax1.set_ylabel('average')
    ax2.set_ylabel('concentration')

    ax1.set_xlabel('phase')
    ax2.legend()

    if SHOW:
        plt.show()
    plt.savefig(file_name + '_phase.pdf', dpi=400)
    plt.close()

    return phase_shift



############# CORRELATION #############
def gauss_detrend(kymo, r):
    global_structs = ndi.gaussian_filter(kymo, sigma=r)
    return kymo - global_structs




def correlate_phase(kymo1, kymo2, file_name, title, upsample_t=1, upsample_x=1,
                    frame_int=1.0, search_range_in_s=50):

    kymo1 = ndi.zoom(kymo1, (upsample_x,upsample_t), order=5)
    kymo2 = ndi.zoom(kymo2, (upsample_x,upsample_t), order=5)

    image_product   = np.fft.fft2(kymo1) * np.fft.fft2(kymo2).conj()
    cc_image        = np.fft.ifft2(image_product)

    unshifted_correlation       = cc_image.real
    time_shifted_correlation    = np.fft.fftshift(unshifted_correlation, axes=1)
    correlation                 = np.fft.fftshift(unshifted_correlation)

    nx, nt = np.shape(correlation)[0], np.shape(correlation)[1]


    ####### Calc MIN #########
    low_b   = int(nt/2 - search_range_in_s/frame_int * upsample_t)
    up_b    = int(nt/2 + search_range_in_s/frame_int * upsample_t)
    min = np.argmin(time_shifted_correlation[0][low_b:up_b]) + low_b

    timestamps = (np.arange(nt) - int(nt/2)) * frame_int/upsample_t
    min_in_s = timestamps[min]


    #### 1D PLOT #####
    correlation1d = time_shifted_correlation[0]

    fig, ax = plt.subplots()

    ax.plot(timestamps, correlation1d/np.max(correlation1d))
    ax.plot(timestamps, kymo1[0], label='radius')
    ax.plot(timestamps, kymo2[0], label='concentration')
    ax.axvline(x=min_in_s, linewidth=1, color='r',
                label='Min: ' + str(np.around(min_in_s, decimals=2)) + ' s')

    ax.set_ylabel('space lag (pixel)')
    ax.set_xlabel('time lag (s)')
    ax.legend()
    plt.savefig(file_name + title + 'correlate1d.pdf', dpi=400)
    if SHOW:
        plt.show()
    plt.close()


    #### IMSHOW PLOT ####
    fig, ax = plt.subplots()

    im = ax.imshow(correlation)
    ax.axvline(x=min, linewidth=0.5, color='r',
                label='Min: ' + str(np.around(min_in_s, decimals=1)) + ' s')

    ax.set_xlabel('time lag (s)')
    ax.set_ylabel('space lag (pixel)')

    ax.set_xticks([int(nt/2)])
    ax.set_xticklabels([0])

    ax.set_yticks([int(nx/2)])
    ax.set_yticklabels([0])

    ax.legend()
    divider = make_axes_locatable(ax)
    fig.colorbar(im, cax = divider.append_axes('right', size='5%', pad=0.05),
                 orientation='vertical')

    fig.tight_layout()

    plt.savefig(file_name + title + 'correlate2d.pdf', dpi=400)
    if SHOW:
        plt.show()
    plt.close()

    return min * frame_int/ upsample_t



def estimate_flow(radii, frame_int):
    radii_b =  bandpass(radii, frame_int, max_freq=0.015)
    integ = np.gradient(radii_b, axis=1) * radii_b

    flow = np.cumsum(integ, axis=0)
    flow /= np.max(flow)

    plt.plot(flow[1])
    plt.plot(flow[-1])

    plt.show()

    show_im(flow[:, 500:-100])

    return


##########################################################################
################################ MAIN ####################################
##########################################################################
SHOW = False
SAVE = False

def main():
    set_keyword     = os.sys.argv[1].strip()
    color           = 'sep'
    method          = 'inter_mean'

    set             = data(set_keyword, no=data(set_keyword).first,
                            method=method, color=color)


    align_keyword   = 'reference_point'

    times           = None, None
    positions       = None, None

    substract_ecto  = False

    max_period      = 140
    range_freq      = 0.001, 0.003
    labels          = get_seeds_positions(set, range_only=True)


    for label in labels:
        kymos_data  = np.load(set.file_dat_set + '_branch_' + str(label) + '.npz')

        path_name = mk_mising_dir(set.file_plot_set + '_branch_' + str(label) + '/')
        file_name = path_name + '/branch_' + str(label) + '_'


        ##################################################################
        ########################## Collect data ##########################
        ##################################################################
        radii               = get_kymo(kymos_data, 'kymo_local_radii',
                                align_keyword, times=times, positions=positions,
                                title='radius(raw)', file_name=file_name, plot=True)

        flow_x              = get_kymo(kymos_data, 'kymo_flow_field_x',
                                align_keyword, times=times, positions=positions,
                                title='flow_x(raw)', file_name=file_name, plot=True)

        flow_y              = get_kymo(kymos_data, 'kymo_flow_field_x',
                                align_keyword, times=times, positions=positions,
                                title='flow_y(raw)', file_name=file_name, plot=True)


        # estimate_flow(radii, set.frame_int)


        kymo_c_green        = get_kymo(kymos_data, 'kymo_c_green',
                                align_keyword, times=times, positions=positions,
                                title='green(raw)', file_name=file_name, plot=True)

        kymo_inner_green    = get_kymo(kymos_data, 'kymo_c_inner_green',
                                align_keyword, times=times, positions=positions)
        kymo_outer_green    = get_kymo(kymos_data, 'kymo_c_outer_green',
                                align_keyword, times=times, positions=positions)

        kymo_c_texas        = get_kymo(kymos_data, 'kymo_c_texas',
                                align_keyword, times=times, positions=positions,
                                title='texas(raw)', file_name=file_name, plot=True)

        kymo_inner_texas    = get_kymo(kymos_data, 'kymo_c_inner_texas',
                                align_keyword, times=times, positions=positions)
        kymo_outer_texas    = get_kymo(kymos_data, 'kymo_c_outer_texas',
                                align_keyword, times=times, positions=positions)


        radii, set.frame_int = shift_radii(radii, set.symm_setup, set.frame_int)
        conce = normalize_green(kymo_c_green, kymo_c_texas, set.symm_setup)
        inner = normalize_green(kymo_inner_green, kymo_inner_texas, set.symm_setup)
        outer = normalize_green(kymo_outer_green, kymo_outer_texas, set.symm_setup)




        plot_kymographs([radii, kymo_c_green, kymo_c_texas, conce, inner, outer],
                        ['radius(cropped)', 'green(cropped)', 'texas(cropped)',
                        'concentration(cropped)', 'inner(cropped)', 'outer(cropped)'],
                        file_name, set.frame_int)


        freq_r = dominant_freq(radii, set.frame_int, max_period=max_period)
        freq_c = dominant_freq(conce, set.frame_int, max_period=max_period)

        print('Dominant frequency: ', freq_r)


        ##################################################################
        #################### Fourier spectrum  ###########################
        ##################################################################
        power_spec(radii, set.frame_int, file_name, 'radius',
                     min_freq=0, max_freq=0.05,
                     mark_freq=freq_r, logscale=True)

        power_spec(radii, set.frame_int, file_name, 'concentration',
                     min_freq=0, max_freq=0.05,
                     mark_freq=freq_r, logscale=True)


        ####################################################
        ########### substract  Ca_global(R) ################
        ####################################################
        if substract_ecto:
            conce = substract_ecto_contribution(radii, conce, set.frame_int)
            inner = substract_ecto_contribution(radii, inner, set.frame_int)
            outer = substract_ecto_contribution(radii, outer, set.frame_int)


        ##################################################################
        ################# calc phase of oscillation ######################
        ##################################################################
        radius_base = bandpass(radii, set.frame_int,
                                min_freq=freq_r-range_freq[0],
                                max_freq=freq_r+range_freq[1])

        conce_base = bandpass(conce, set.frame_int,
                                min_freq=freq_r-range_freq[0],
                                max_freq=freq_r+range_freq[1])

        phase_radius, _, _ = extract_phase(radius_base, file_name, 'radius', set.frame_int)
        phase_conce, _, _  = extract_phase(conce_base, file_name, 'concentration', set.frame_int)



        ##################################################################
        ############ bandpass kymographs before binning ##################
        ##################################################################
        radii_b = bandpass(radii, set.frame_int,  min_freq=freq_r-range_freq[0],
                                                max_freq=freq_r+range_freq[1])
        conce_b = bandpass(conce, set.frame_int,  min_freq=freq_r-range_freq[0],
                                                max_freq=freq_r+range_freq[1])
        inner_b = bandpass(inner, set.frame_int,  min_freq=freq_r-range_freq[0],
                                                max_freq=freq_r+range_freq[1])
        outer_b = bandpass(outer, set.frame_int,  min_freq=freq_r-range_freq[0],
                                                max_freq=freq_r+range_freq[1])


        ##################################################################
        ######################## Phase dependence ########################
        ##################################################################
        titles = ['radius(band)', 'concentration(band)', 'inner(band)', 'outer(band)']
        colors = ['orange', 'blue', 'darkblue', 'lightskyblue']
        kymos  = [radii_b, conce_b, inner_b, outer_b]

        plot_kymographs(kymos, titles, file_name, set.frame_int)


        phase_shift = phase_average(phase_radius, kymos, titles, colors,
                                    file_name)


        ##################################################################
        ######################### correlations ###########################
        ##################################################################
        min_k = correlate_phase(radii_b, conce_b, file_name, 'kymos',
                                upsample_t=10, upsample_x=1,
                                frame_int=set.frame_int)


        min_p = correlate_phase(phase_radius, phase_conce, file_name, 'phases',
                                upsample_t=10, upsample_x=1,
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
                    out_var.write('# data_set \t label \t' +
                                    'phase_radius \t phase_conc \t' +
                                    'phase_inner \t phase_outer \t' +
                                    'min_kymo + \t min_phase \n')

            with open(data_sets_summary, "a") as out_var:
                out_var.write(set_keyword +
                                '\t' + str(label) +
                                '\t' + str(phase_shift[0]) +
                                '\t' + str(phase_shift[1]) +
                                '\t' + str(phase_shift[2]) +
                                '\t' + str(phase_shift[3]) +
                                '\t' + str(min_k) +
                                '\t' + str(min_p) + '\n')


if __name__ == '__main__':
    main()
