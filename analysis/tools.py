import numpy as np
import datetime

def log_message(set, script, *add_args):
    log_message = str(datetime.datetime.now()) + ': ' + script
    for arg in add_args:
        log_message += ' ' + str(arg)
    print(log_message, file=open(set.file_log, "a"))


def closest_number(arr, ref):
    return arr[np.argmin((arr - ref)**2)]


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
