import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.signal import hilbert
from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import curve_fit
from scipy import signal


import os
import itertools
import copy


from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *


from skimage import feature
from skimage import transform
from skimage.util.shape import view_as_windows
from skimage.feature import match_template


class kymograph:
    """
    class that collects kymograph data array, axis information and useful
    information for plotting:

        self.kymo       data array

        self.name       name of the quantity eg radius,
        self.state      denotes current state of the analysis

        self.t_scaling  scalings used for axis labels
        self.x_scaling  ...
        self.c_scaling  ...

        self.cmap       cmap used for plotting
        self.color      color used for plotting


    """

    # ===================================================== #
    #  Constructor (-like) and copy methods  #
    # ===================================================== #
    def __init__(self, kymo, name, t_scaling, x_scaling, c_scaling,
                state='raw', cmap='viridis', color=None):

        self.kymo       = kymo

        self.name       = name
        self.state      = state

        self.t_scaling  = t_scaling
        self.x_scaling  = x_scaling
        self.c_scaling  = c_scaling

        self.cmap       = cmap
        self.color      = color



    def copy(self):
        return copy.deepcopy(self)


    def copy_meta(self, new_kymo):
        return kymograph(new_kymo, self.name,
                            self.t_scaling,
                            self.x_scaling,
                            self.c_scaling,
                            self.state,
                            self.cmap,
                            self.color)


    def get_title(self, for_file_name=True):
        if for_file_name:
            return self.name.replace(' ', '') + '(' + self.state.replace(' ', '') + ')'

        else:
            if self.state == 'raw':
                return self.name
            elif self.state == 'cropped':
                return 'cropped ' + self.name
            elif self.state == 'band':
                return 'bandpassed ' + self.name


    # ===================================================== #
    # Operations on existing data, returning new instances
    # of kymograph with state 'raw', 'cropped' or 'band'
    # ===================================================== #

    def get_dat(kymos_data, keyword,
                name, t_s, x_s, c_s, cmap,
                align_keyword='reference_point',
                times=(None,None), positions=(None, None)):

        alignment = kymos_data['alignment']

        kymo      = kymos_data[keyword]

        if align_keyword == None:
            kymo = np.transpose(kymo)
        else:
            kymo = align_kymo(kymo, align_keyword, alignment=alignment)
            kymo = np.transpose(kymo)

        kymo = kymo[slice(*positions), slice(*times)]

        return kymograph(kymo, name, t_s, x_s, c_s, state='raw', cmap=cmap, color=None)


    def crop_kymo(self):
        crop = crop_nans(self.kymo)

        new = self.copy_meta(crop)
        new.state = 'cropped'
        return new


    def bandpass(self, min_freq=None, max_freq=None,
                min_period=None, max_period=None):

        # use period if given:
        if max_period != None:
            min_freq = 1./max_period

        if min_period != None:
            max_freq = 1./min_period

        kymo_f = np.fft.rfft(self.kymo)
        freq    = np.fft.rfftfreq(self.kymo.shape[-1], d=self.t_scaling)

        if max_freq != None:
            kymo_f[:, (freq > max_freq)] = 0
        if min_freq != None:
            kymo_f[:, (freq < min_freq)] = 0

        band = np.fft.irfft(kymo_f)
        new = self.copy_meta(band)
        new.state = 'band'
        return new



    # ===================================================== #
    # Create kymograph that is used in signal analysis,
    # includes cropping and shifting if setup is symm(etric) #
    # ===================================================== #
    def shift(self, symm_setup, color):
        new = self.crop_kymo()
        if symm_setup:
            shifted = ndi.zoom(new.kymo, (1,2), order=5)[1:]
            new = new.copy_meta(shifted)
            new.t_scaling = new.scaling/2.

        else:
            new = self.crop_kymo()
            new.color = color
        return new


    def ca_concentration(green, texas, symm_setup, name, cmap, color):
        green_c = green.crop_kymo()
        texas_c = texas.crop_kymo()

        if symm_setup:
            g2 = ndi.zoom(green_c.kymo, (1,2), order=5)[:-1]
            t2 = ndi.zoom(texas_c.kymo, (1,2), order=5)[1:]

            new = kymograph(g2/t2, name,
                            green_c.t_scaling/2.,
                            green_c.x_scaling,
                            green_c.c_scaling,
                            green_c.state,
                            cmap)

        else:
            new = kymograph(green_c.kymo/texas_c.kymo, name,
                            green_c.t_scaling,
                            green_c.x_scaling,
                            green_c.c_scaling,
                            green_c.state,
                            cmap)
        new.state = 'cropped'
        new.color = color

        return new





# ===================================================== #
#  Plotting  #
# ===================================================== #

def plot_kymograph(kymo_list, filename, lim_mode='minmax'):

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

    for kymo in kymo_list:
        fig, ax = plt.subplots()
        if lim_mode != 'minmax':
            std_im, mean_im = np.nanstd(kymo.kymo), np.nanmean(kymo.kymo)
            min, max = mean_im - lim_mode*std_im, mean_im + lim_mode*std_im

            im = ax.imshow(kymo.kymo, cmap=kymo.cmap, aspect='auto',
                            vmin = min, vmax = max)
        else:
            im = ax.imshow(kymo.kymo, cmap=kymo.cmap, aspect='auto')

        ax.set_xlim(left=0, right=np.shape(kymo.kymo)[1])
        locs = ax.get_xticks()
        ax.set_xticks(locs)
        ax.set_xticklabels(np.around(locs * kymo.t_scaling, decimals=0).astype(int))

        ax.set_title(kymo.get_title(0))
        ax.set_ylabel('space (pixel)')
        ax.set_xlabel('time (s)')

        divider = make_axes_locatable(ax)
        fig.colorbar(im, cax = divider.append_axes('right', size='5%', pad=0.05),
                     orientation='vertical')

        plt.savefig(filename + kymo.get_title() + '.pdf', dpi=400)
        plt.close()




# ===================================================== #
#  time series  #
# ===================================================== #

def slid_aver(arr, window_size):
    pad_size = window_size
    expanded = np.pad(arr, pad_size , mode='reflect')
    new = np.convolve(expanded, np.ones(window_size)/window_size, mode='same')
    return new[pad_size:-pad_size]


def slid_aver2d(arr2d, window_size, axis=0):
    return np.apply_along_axis(slid_aver, axis, arr2d, *(window_size,))


def plot_time_series(kymos, filename=None, window_size=10):

    fig, axes = plt.subplots(len(kymos),1, sharex=False)
    frame_int = kymos[0].t_scaling

    for ax, kymo in zip(axes, kymos):

        time_series = slid_aver2d(kymo.kymo, window_size)[::window_size]
        timestamps = np.arange(np.shape(time_series)[1]) * frame_int


        pixel_positions = np.arange(0, np.shape(time_series)[0] ) * window_size

        norm = mpl.colors.Normalize(vmin=pixel_positions.min(),
                                    vmax=pixel_positions.max())

        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='plasma')
        cmap.set_array([])

        for t_s, p in zip(time_series, pixel_positions):
            im = ax.plot(timestamps, t_s, c=cmap.to_rgba(p))

        ax.set_xlim(0, timestamps[-1])
        ax.set_title(kymo.get_title(0))

        if not ax == axes[-1]:
            ax.set_xticks([])


    plt.savefig(filename + 'time_series.pdf', dpi=400, bbox_inches='tight')
    plt.close()
    # plt.show()


# ===================================================== #
#  FFT  #
# ===================================================== #


def dominant_freq(kymo, min_freq=0.05, max_period=None):

    if max_period!=None:
        min_freq = 1./max_period

    f, Pxx_den_full = signal.periodogram(kymo.kymo, fs=1/kymo.t_scaling)
    Pxx_den = np.mean(Pxx_den_full, axis=0)

    Pxx_den = Pxx_den[f>min_freq]
    f = f[f>min_freq]


    dominant_freq = f[np.argmax(Pxx_den)]
    return dominant_freq



def power_spec(kymo, filename,
                min_freq=0.005, max_freq=0.05,          #periods between 200 and 20s
                min_period=None, max_period=None,
                mark_freq=None, logscale=False):


    if min_period!=None and max_period!=None:
        min_freq = 1./max_period
        max_freq = 1./min_period

    f1, Pxx_den1 = signal.periodogram(kymo.kymo, fs=1./kymo.t_scaling,
                                scaling='spectrum')
    # Pxx_den1 = ndi.gaussian_filter(Pxx_den1, sigma=2)


    range = (f1 < max_freq) * (f1 > min_freq)

    Pxx_den1    = Pxx_den1[:,range]
    f1          = f1[range] * 1e3 ### go to mHz

    if logscale:
        Pxx_den1 = np.log(Pxx_den1)

    fig, axes = plt.subplots(2,1)
    ax = axes.ravel()



    im = ax[0].imshow(Pxx_den1, aspect='auto', cmap='plasma')

    cax = fig.add_axes([0.2, 0.95, 0.6, 0.05])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    if logscale:
        cbar.ax.set_title('log(power)', fontsize='medium')
    else:
        cbar.ax.set_title('power', fontsize='medium')

    # ax[0].set_title('log(power) ' + title)
    ax[0].set_xticks([])
    ax[0].set_ylabel('pixel')


    # divider = make_axes_locatable(ax[0])
    # fig.colorbar(im, cax = divider.append_axes('right',size='5%', pad=0.05),
    #              orientation='horizontal')
    # fig.colorbar(im, cax = ax[0], orientation='horizontal', fraction=.1)



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
    # fig.tight_layout()


    plt.savefig(filename + kymo.get_title() + '_power.pdf', dpi=200, bbox_inches='tight')
    plt.close()



# ===================================================== #
#  Phase dependence  #
# ===================================================== #

#  hilbert  #

def extract_phase(kymo, filename, min_freq=None, max_freq=None):

    if min_freq!=None and max_freq!=None:
        kymo_b = kymo.bandpass(min_freq, max_freq)
    else:
        kymo_b = kymo.copy()


    analytic_signal = hilbert(kymo_b.kymo)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.angle(analytic_signal)
    # instantaneous_phase = np.unwrap(instantaneous_phase, axis=1)
    instantaneous_frequency = np.diff(instantaneous_phase)/(2.0*np.pi)*kymo_b.t_scaling

    fig, axes = plt.subplots(1,3, figsize=(6, 3), sharex=True, sharey=True)
    ax = axes.ravel()
    ax[0].set_title(kymo_b.get_title(0))
    c = ax[0].imshow(kymo_b.kymo, aspect='auto')

    ax[1].set_title('amplitude')
    ax[1].imshow(amplitude_envelope, cmap='plasma', aspect='auto')

    ax[2].set_title('phase')
    ax[2].imshow(instantaneous_phase, cmap='twilight_shifted', aspect='auto')

    # ax[3].set_title('frequency')
    # ax[3].imshow(instantaneous_frequency)

    plt.savefig(filename + kymo_b.get_title() + 'hilbert.pdf', dpi=400)


    plt.close()

    return instantaneous_phase, amplitude_envelope, instantaneous_frequency


# binning #
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


def phase_average(phase, kymos, filename, y_label):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    phase_shift = []

    for k in kymos:
        bin_mean, bin_edges, bin_std = bin_acc_phase(phase, k.kymo, n_bins=15,
                                                    norm=False)
        bin = bin_edges[:-1] + np.diff(bin_edges)/2


        phase_sample = np.linspace(bin[0], bin[-1], 100)
        popt, pcov   = curve_fit(sin_func, bin, bin_mean, p0=(-1,0))

        if k == kymos[0]:
            axis = ax1
        else:
            axis = ax2

        axis.plot(bin, bin_mean , label=k.get_title(0) + ', shift: ' +
                    str( np.around(popt[-1], decimals=3) ), color=k.color)

        # plt.fill_between(bin, bin_mean - bin_std, bin_mean + bin_std, alpha=0.2, color=k.color)

        # min = bin[np.argmin(bin_mean)]
        # plt.axvline(x=min, linewidth=1, color=k.color)

        axis.plot(phase_sample, sin_func(phase_sample, *popt), '--', color=k.color)
        y1 = sin_func(phase_sample, popt[0], popt[1] + pcov[1,1]**0.5)
        y2 = sin_func(phase_sample, popt[0], popt[1] - pcov[1,1]**0.5)
        # plt.fill_between(phase_sample, y1, y2, color=k.color, alpha=0.15)

        axis.axvline(x=popt[-1], linewidth=1, color=k.color)
        phase_shift.append(popt[-1])


    ax1.set_ylabel('radius')
    ax2.set_ylabel(y_label)

    ax1.set_xlabel('phase')
    ax2.legend()


    plt.savefig(filename +'_'+ y_label + '_phase.pdf', dpi=400, bbox_inches='tight')
    plt.close()

    return phase_shift




# ===================================================== #
#  Correlation  #
# ===================================================== #

def gauss_detrend(kymo, r):
    global_structs = ndi.gaussian_filter(kymo, sigma=r)
    return kymo - global_structs



def correlate_phase(kymo_a, kymo_b, filename, title, upsample_t=1, upsample_x=1,
                    search_range_in_s=50, detrend=None):


    kymo1 = kymo_a.copy()
    kymo2 = kymo_b.copy()

    if detrend=='gauss':
        sec = 150
        sigma_gauss = sec / kymo1.t_scaling
        kymo1.kymo = gauss_detrend(kymo1.kymo, sigma_gauss)
        kymo2.kymo = gauss_detrend(kymo2.kymo, sigma_gauss)


    kymo1.kymo = ndi.zoom(kymo1.kymo, (upsample_x,upsample_t), order=5)
    kymo2.kymo = ndi.zoom(kymo2.kymo, (upsample_x,upsample_t), order=5)

    image_product   = np.fft.fft2(kymo1.kymo) * np.fft.fft2(kymo2.kymo).conj()
    cc_image        = np.fft.ifft2(image_product)

    unshifted_correlation       = cc_image.real
    time_shifted_correlation    = np.fft.fftshift(unshifted_correlation, axes=1)
    correlation                 = np.fft.fftshift(unshifted_correlation)

    nx, nt = np.shape(correlation)[0], np.shape(correlation)[1]


    ####### Calc MIN #########
    low_b   = int(nt/2 - search_range_in_s/kymo1.t_scaling * upsample_t)
    up_b    = int(nt/2 + search_range_in_s/kymo1.t_scaling * upsample_t)
    min = np.argmin(time_shifted_correlation[0][low_b:up_b]) + low_b

    timestamps = (np.arange(nt) - int(nt/2)) * kymo1.t_scaling/upsample_t
    min_in_s = timestamps[min]


    #### 1D PLOT #####
    correlation1d = time_shifted_correlation[0]

    fig, ax = plt.subplots()

    ax.plot(timestamps, correlation1d/np.max(correlation1d))
    ax.plot(timestamps, kymo1.kymo[0], label='radius')
    ax.plot(timestamps, kymo2.kymo[0], label='concentration')
    ax.axvline(x=min_in_s, linewidth=1, color='r',
                label='Min: ' + str(np.around(min_in_s, decimals=2)) + ' s')

    ax.set_ylabel('space lag (pixel)')
    ax.set_xlabel('time lag (s)')
    ax.legend()
    plt.savefig(filename + title + 'correlate1d.pdf', dpi=400)

    plt.close()


    #### IMSHOW PLOT ####
    fig, ax = plt.subplots()

    im = ax.imshow(correlation, aspect='auto')
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

    plt.savefig(filename + title + 'correlate2d.pdf', dpi=400)

    plt.close()

    return min * kymo1.t_scaling/ upsample_t



# ===================================================== #
#  Time/space dependend phase_shift #
# ===================================================== #
def phase_shift2d(kymo_a, kymo_b, filename, upsample=1, window_size=50,
                        sampling=20, search_range_in_s=40):

    """
    upsample            upsample images to get subpixel precision
    window_size         window for corr given by window_size * window_size
    sampling            sampling step between pixels that serve as seeds for window
    search_extend       extends search area by search_extend
    """

    kymo1 = kymo_a.copy()
    kymo2 = kymo_b.copy()

    search_extend = np.around(search_range_in_s/kymo1.t_scaling).astype(int)

    s_min = search_extend
    s_max = -search_extend

    window  = (window_size + 2*search_extend * upsample,)*2 # calc window shape eg (40,40)
    step    = sampling * upsample


    if upsample > 1:
        kymo1.kymo = ndi.zoom(kymo1.kymo, upsample, order=5)
        kymo2.kymo = ndi.zoom(kymo2.kymo, upsample, order=5)


    a_windows = view_as_windows(kymo1.kymo, window_shape= window, step=step)
    b_windows = view_as_windows(kymo2.kymo, window_shape= window, step=step)
    shape = np.shape(b_windows)

    mins = np.zeros(shape[:2])

    for i in range(shape[0]):
        for j in range(shape[1]):
            # crop window to orig. window_size
            a = a_windows[i,j][s_min:s_max, s_min:s_max]
            b = b_windows[i,j]

            corr2d = match_template(b,a)
            corr1d = corr2d[int(corr2d.shape[0]/2)]

            plt.plot(corr1d)
            plt.show()

            ind = np.unravel_index(np.argmin(corr1d, axis=None), corr1d.shape)
            mins[i,j] = ind[0] - int(corr1d.shape[0]/2)

    fig, ax = plt.subplots()
    ax.imshow(mins, aspect='auto')
    plt.savefig(filename  + 'phase_shift_map.pdf', dpi=400)

    return mins * kymo1.t_scaling/ upsample



def estimate_flow(radii, frame_int):
    radii_b =  bandpass(radii, frame_int, max_freq=0.015)
    integ = np.gradient(radii_b, axis=1) * radii_b

    flow = np.cumsum(integ, axis=0)
    flow /= np.max(flow)

    plt.plot(flow[1])
    plt.plot(flow[-1])

    # plt.show()

    show_im(flow[:, 500:-100])

    return
