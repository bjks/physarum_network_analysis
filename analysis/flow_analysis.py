import numpy as np
import matplotlib.pyplot as plt

from skimage import feature
from skimage import transform
from skimage.util.shape import view_as_windows

from scipy import ndimage as ndi
from scipy import signal as signal


def signed_avg(a, b, for_sign=1):
    return np.sqrt( np.dot(a,a) + np.dot(b,b) ) * np.sign(for_sign)



def bring_to_shape_of(ar, ar_with_desired_shape, order=0):
    """ expands array to desired shape"""
    return transform.resize(ar, np.shape(ar_with_desired_shape),
                            order=order, mode='reflect')



def flow_quantification(frame_a, frame_b, mask, upsample=1, window_size=50,
                        sampling=20, search_extend=20, return_scalar_v=True):
    """
    mask                windows that contain only masked pixel are ignored
    upsample            upsample images to get subpixel precision
    window_size         window for corr given by window_size * window_size
    sampling            sampling step between pixels that serve as seeds for window
    search_extend       extends search area by search_extend
    return_scalar_v     returns (signed) scalar value for velocity
    """
    if search_extend<=0:
        print('WARNING: search_extend <=0! Flow will be zero!')

    window  = (window_size + search_extend * upsample,)*2 # calc window shape eg (40,40)
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
            a = a_windows[i,j][search_extend:-search_extend, search_extend:-search_extend]
            m = m_windows[i,j][search_extend:-search_extend, search_extend:-search_extend]


            # must contain non masked pixels. must not be zeros only to avoid artifacts
            if np.sum(m)>0 and np.sum(a)>0:
                b = b_windows[i,j]
                corr = signal.correlate2d(a, b, mode='valid')

                corr_shape = corr.shape
                shift = np.unravel_index(np.argmax(corr, axis=None), corr_shape)

                flow_field_x[i,j] = shift[1]  - int(corr_shape[1]/2)
                flow_field_y[i,j] = shift[0]  - int(corr_shape[0]/2)



    flow_field_x/=upsample  # rescale flow to orignal pixel distance
    flow_field_y/=upsample

    # plt.quiver(flow_field_x, flow_field_y)
    # plt.show()


    if return_scalar_v:
        flow_x = np.nanmean(flow_field_x)
        flow_y = np.nanmean(flow_field_y)
        avg_flow = signed_avg(flow_x, flow_y, for_sign = flow_y)
        return avg_flow, flow_field_x, flow_field_y

    return flow_field_x, flow_field_y



def average_flow_over_frames(frames, mask, upsample=1, window_size=50,
                                sampling=20, search_extend=5):

    flow_field_x = []
    flow_field_y = []
    for j in range(len(frames)-1):
        x, y, = flow_quantification(frames[j], frames[j+1], mask, upsample=upsample,
                                window_size=window_size, sampling=sampling,
                                search_extend=search_extend,
                                return_scalar_v=False)


        flow_field_x.append(x)
        flow_field_y.append(y)

    mean_flow_x = np.mean(flow_field_x, axis=0)
    mean_flow_y = np.mean(flow_field_y, axis=0)

    plt.quiver(mean_flow_x, mean_flow_y)
    plt.show()

    return  bring_to_shape_of(mean_flow_x, frames[0]), \
            bring_to_shape_of(mean_flow_y, frames[0])
