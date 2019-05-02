import numpy as np


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


def crop_aligned_kymo(kymo):
    new_kymo = []
    for arr in np.transpose(kymo)[:]:
        if np.logical_not(np.isnan(arr).any()):
            new_kymo.append(arr)

    return np.transpose(new_kymo)
