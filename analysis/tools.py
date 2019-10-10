import numpy as np

def fft_kymo(signal, frame_int):

    kymo_f  = np.fft.rfft(signal)

    freq    = np.fft.rfftfreq(signal.shape[-1]) / frame_int
    fourier = np.mean(kymo_f.real, axis=0)
    plt.plot(freq, fourier)
    plt.show()


def fill_up_array(array, add1, add2):

    new_array = np.append(np.full(add1, np.nan), array)
    new_array = np.append(new_array, np.full(add2, np.nan))

    return new_array


def align_kymo(kymo, mode='start', alignment=None):
    new_kymo = []
    if mode =='start':
        max_len = np.max([len(arr) for arr in kymo[:]])
        for arr in kymo[:]:
            add = max_len - len(arr)
            new_kymo.append(fill_up_array(arr, 0, add))

    elif mode =='final':
        max_len = np.max([len(arr) for arr in kymo[:]])
        for arr in kymo[:]:
            add = max_len - len(arr)
            new_kymo.append(fill_up_array(arr, add, 0))

    elif mode =='reference_point':
        max_al  = np.max(alignment)
        temp = []
        for arr, al in zip(kymo[:], alignment):
            add1 = max_al - al
            temp.append(fill_up_array(arr, add1, 0))

        max_len = np.max([len(arr) for arr in temp[:]])
        for arr in temp[:]:
            add = max_len - len(arr)
            new_kymo.append(fill_up_array(arr, 0, add))

    elif mode == 'center':
        max_len = np.max([len(arr) for arr in kymo[:]])
        for arr in kymo[:]:
            add = max_len - len(arr)
            add1 = (add/2).astype(int)
            add2 = add - add1
            new_kymo.append(fill_up_array(arr, add1, add2))

    return np.array(new_kymo)


def crop_nans(arr2d):
    new_arr2d = []
    for arr in arr2d[:]:
        if np.logical_not(np.isnan(arr).any()):
            new_arr2d.append(arr)

    return np.array(new_arr2d)
