import numpy as np
import matplotlib.pyplot as plt

from skimage import feature
from skimage import transform
from skimage.util.shape import view_as_windows
from skimage.feature import match_template
# from analysis.network_analysis import *

from scipy import ndimage as ndi
from scipy import signal as signal
import time


def show_flow(x,y):
    plt.quiver(x,y, scale=1, scale_units='x')
    plt.show()


def signed_avg(a, b, for_sign=1):
    return np.sqrt( np.dot(a,a) + np.dot(b,b) ) * np.sign(for_sign)


def bring_to_shape_of(ar, ar_with_desired_shape, order=0):
    """ expands array to desired shape"""
    return transform.resize(ar, np.shape(ar_with_desired_shape),
                            order=order, mode='wrap', anti_aliasing=True)


def mask_filter(arr, mask=None, size=3, selem='square'):
    if np.size(mask)>1:
        arr_filtered = np.where(mask!=0, arr, np.nan) # 0 -> NaN (to be ignored)
    else:
        arr_filtered = arr

    if selem == 'square':
        kernel = np.ones((size,size))
    elif 'disk':
        kernel = morph.disk(size)

    arr_filtered = convolve_fft(arr_filtered, kernel, fill_value=np.nan)
    return np.nan_to_num(arr_filtered) # NaN -> 0


def replace_outlier_fq(arr, thresh, size=3, mask=None, replace_with='mean'):

    if np.size(mask)>1: # ... to allow mask = None
        background = mask_filter(arr, mask=mask, size=size)
        outlier_mask   = np.where(arr < thresh * background, 1., 0.) * mask
        estim_back   = mask_filter(arr, mask=outlier_mask, size=size)
        corrected_arr = np.where(outlier_mask!=0, arr, estim_back) * mask

    else:
        background   = mask_filter(arr, size=size)
        outlier_mask   = np.where(arr < thresh * background, 1., 0.)
        estim_back   = mask_filter(arr, mask=outlier_mask, size=size)
        corrected_arr = np.where(outlier_mask!=0, arr, estim_back)

    return corrected_arr



def flow_quantification(frame_a, frame_b, mask,
                        upsample=1,
                        window_size=50,
                        sampling=20,
                        search_extend=20,
                        return_scalar_v=True,
                        corr_tresh=0.5,
                        outlier_thresh=None):
    """
    mask                windows that contain only masked pixel are ignored
    upsample            upsample images to get subpixel precision
    window_size         window for corr given by window_size * window_size
    sampling            sampling step between pixels that serve as seeds for window
    search_extend       extends search area by search_extend
    return_scalar_v     returns (signed) scalar value for velocity
    corr_tresh          minimal max(corr) to be accepted for flow field
    outlier_thresh      replace outlier with local mean if value > outlier_thresh * local mean
    """

    s_min = search_extend * upsample
    s_max = - search_extend * upsample

    window  = ( (window_size + 2*search_extend) * upsample ,)*2 # calc window shape eg (40,40)
    step    = sampling * upsample


    if upsample > 1:
        mask    = ndi.zoom(mask.astype(bool), upsample, order=0) # order 0, to keep bool behavior!!
        frame_a = ndi.zoom(frame_a, upsample, order=5)
        frame_b = ndi.zoom(frame_b, upsample, order=5)


    a_windows = view_as_windows(frame_a, window_shape= window, step=step)
    b_windows = view_as_windows(frame_b, window_shape= window, step=step)
    m_windows = view_as_windows(mask, window_shape= window, step=step)
    shape = np.shape(b_windows)

    flow_field_x = np.zeros(shape[:2])
    flow_field_y = np.zeros(shape[:2])

    for i in range(shape[0]):
        for j in range(shape[1]):
            # crop window to orig. window_size
            a = a_windows[i,j][s_min:s_max, s_min:s_max]
            m = m_windows[i,j][s_min:s_max, s_min:s_max]


            # must only contain non masked pixels, alternatively use .any()
            if np.all(m):
                b = b_windows[i,j]

                corr = match_template(b,a)

                if np.max(corr) > corr_tresh and np.any(a):

                    corr_shape = corr.shape
                    shift = np.unravel_index(np.argmax(corr, axis=None), corr_shape)

                    flow_field_x[i,j] = shift[1] - int(corr_shape[1]/2)
                    flow_field_y[i,j] = shift[0] - int(corr_shape[0]/2)


    flow_field_x /= upsample
    flow_field_y /= upsample

    if outlier_thresh != None:
        flow_field_x = replace_outlier_fq(flow_field_x, outlier_thresh)
        flow_field_y = replace_outlier_fq(flow_field_y, outlier_thresh)


    if return_scalar_v:
        flow_x = np.nanmean(flow_field_x)
        flow_y = np.nanmean(flow_field_y)
        avg_flow = signed_avg(flow_x, flow_y, for_sign = flow_y)
        return avg_flow, flow_field_x, flow_field_y

    return flow_field_x, flow_field_y



def average_flow_over_frames(frames, mask,
                            upsample=1,
                            window_size=50,
                            sampling=20,
                            search_extend=20,
                            corr_tresh=0.5,
                            outlier_thresh=None):
    """
    averages over [ flow_quantification(frame[i], frame[i+1]),
                    flow_quantification(frame[i+1], frame[i+2]) ... ]


    tunnels params to flow_quantification

    returns flow_field_x/y with the same shape as frames[0]

    """
    flow_field_x = []
    flow_field_y = []
    for j in range(len(frames)-1):

        x, y, = flow_quantification(frames[j], frames[j+1], mask,
                                    upsample=upsample,
                                    window_size=window_size,
                                    sampling=sampling,
                                    search_extend=search_extend,
                                    return_scalar_v=False,
                                    corr_tresh=corr_tresh,
                                    outlier_thresh=outlier_thresh)

        flow_field_x.append(x)
        flow_field_y.append(y)



    mean_flow_x = np.mean(flow_field_x, axis=0)
    mean_flow_y = np.mean(flow_field_y, axis=0)
    # show_flow(mean_flow_x, mean_flow_y)


    return  bring_to_shape_of(mean_flow_x, frames[0]), \
            bring_to_shape_of(mean_flow_y, frames[0])
