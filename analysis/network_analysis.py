import numpy as np
import matplotlib.pyplot as plt

import sys

from skimage.feature import peak_local_max
import skimage.morphology as morph
import skimage.filters.rank as rank
from skimage.measure import regionprops


from scipy import ndimage as ndi
from scipy import LowLevelCallable
from scipy.ndimage.filters import generic_filter


import numba
from numba import cfunc, carray
from numba.types import intc, CPointer, float64, intp, voidptr

############ File reading ########
def greyscale(image):
    return np.dot(image[...,:3], [0.299, 0.587, 0.114])

def read_file(filename):
    image = plt.imread(filename, 'tif').astype(float)
    if image.ndim==3:
        return greyscale(image)
    else:
        return image

### create diluted skeleton and masks the background for plotting ###
def thick_skeleton(skeleton, times = 10):
    selem = morph.disk(times)
    thick_skeleton = skeleton.copy()

    thick_skeleton = morph.dilation(skeleton, selem)
    thick_skeleton = np.ma.masked_where(thick_skeleton == 0, thick_skeleton)
    return thick_skeleton

#####################################
############# Creating mask #########
#####################################

def create_mask(dye, sig, thresh, halo_sig=None):
    ### if halo_sig is given the local background is estimated and used for mask
    if halo_sig != None:
        halo = ndi.gaussian_filter(dye, sigma=halo_sig)
        mask = ndi.gaussian_filter(dye, sigma=sig)
        mask = np.where(mask > thresh * halo, 1., 0.)
    ### pixels are compared to mean of entire image
    else:
        mask = ndi.gaussian_filter(dye, sigma=sig)
        mask = np.where(mask > thresh*(np.mean(dye)), 1., 0.)
    return mask

def extract_nerwork(mask):
    labels = morph.label(mask, connectivity=2)
    regions = regionprops(labels)

    label_max = regions[np.argmax([r.area for r in regions])].label
    return np.where(labels == label_max, 1., 0.)

#####################################
######## image interpolation ########
#####################################
def interpl_masks(mask1, mask2):
    mask1 = mask1.astype(bool)
    mask2 = mask2.astype(bool)

    d1 = -ndi.distance_transform_edt(mask1) + ndi.distance_transform_edt(np.invert(mask1))
    d2 = -ndi.distance_transform_edt(mask2) + ndi.distance_transform_edt(np.invert(mask2))

    d = np.mean([d1, d2], axis=0)
    new_mask = np.where(d < 0, 1, 0)
    return new_mask

def interpl_dye(raw1, raw2):
    return np.mean([raw1, raw2], axis = 0)

def interpl_images(dye1, dye2, sigma, threshold, halo_sig):

    mask1 = create_mask(dye1, sigma, threshold, halo_sig)
    mask2 = create_mask(dye2, sigma, threshold, halo_sig)

    mask = interpl_masks(mask1, mask2)
    interpl = interpl_dye(dye1, dye2)
    return interpl * mask, mask

#####################################
########## skeleton mapping #########
#####################################

def extract_skeleton(mask):
    # temp = ndi.gaussian_filter(dye, sigma=sig)
    # temp = np.where(temp > thresh*(np.mean(temp)) , 1., 0.)
    return morph.skeletonize(mask)

def extract_radii(mask, skeleton):
    distance = ndi.distance_transform_edt(mask)
    local_radii = distance * skeleton
    return local_radii

#########################################################################
########## replace spots with estim. background, based on nearby pixels #
#                           jit stuff is just for performance   #########
#########################################################################

# speeding up generic filter
def jit_filter_function(filter_function):
    jitted_function = numba.jit(filter_function, nopython=True)
    @cfunc(intc(CPointer(float64), intp, CPointer(float64), voidptr))
    def wrapped(values_ptr, len_values, result, data):
        values = carray(values_ptr, (len_values,), dtype=float64)
        result[0] = jitted_function(values)
        return 1
    return LowLevelCallable(wrapped.ctypes)

@jit_filter_function
def nan_mean(values):
    return np.nanmean(values)
###

def disk_filter(dye, mask, r):
    dye_filtered = np.where(mask!=0, dye, np.nan) # 0 -> NaN (to be ignored in nanmean)
    selem        = morph.disk(r)

    dye_filtered = generic_filter(dye_filtered, nan_mean, footprint=selem,
                                  mode='constant', cval=np.nan)

    return np.nan_to_num(dye_filtered) # NaN -> 0

def remove_spots(dye, mask, spots_sig, thresh_spots):
    #### calculate background ####
    selem        = morph.disk(spots_sig)
    background   = disk_filter(dye, mask, spots_sig)

    #### calculate mask ####
    spots_mask   = np.where(dye < thresh_spots * background, 1., 0.) * mask

    #### replace spots with estimated background ###
    estim_back   = disk_filter(dye, spots_mask, spots_sig)
    dye_spotless = np.where(spots_mask!=0, dye, estim_back ) * mask

    return spots_mask, dye_spotless



#####################################
########## intensity calcs #########
#####################################

def calc_ratio(dye, ref):
    return np.true_divide(dye, ref, out=np.zeros_like(dye), where=ref!=0)

#####################################
########## skeleton mapping #########
#####################################

def tube_radius_at_point(mask, skel, local_radii):
    ### calculate the radius (or other given quantity) at every position
    #   in network based on value of the nearest skeleton point
    radii = mask.copy().astype(float)
    inds = ndi.distance_transform_edt(np.invert(skel.astype(bool)),
                                      return_distances=False, return_indices=True)

    radii = radii.flatten()
    inds0 = inds[0].flatten()
    inds1 = inds[1].flatten()

    for i in  range(0, len(radii)):
        if radii[i]!=0:
            radii[i] = local_radii[inds0[i]][inds1[i]]

    return radii.reshape(np.shape(mask))

def relative_distance(skeleton, mask, local_radii):
    distance = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                          return_distances=True, return_indices=False)
    distance *= mask
    radii = tube_radius_at_point(mask, skeleton, local_radii)
    return np.true_divide(distance, radii, out=np.zeros_like(radii), where=radii!=0)

def absolute_distance(skeleton, mask):
    distance = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                          return_distances=True, return_indices=False)
    return distance


def circle_mean(dye, skeleton, mask, local_radii):
    cx = np.arange(0, dye.shape[1])
    cy = np.arange(0, dye.shape[0])

    #### calc realtive distance (r/R) to skeleton at every position
    relative_dist = relative_distance(skeleton, mask, local_radii)

    ### seperate outer and inner of network based on relative_dist
    dye_inner = np.where(relative_dist <= 0.6, dye, 0)
    dye_outer = np.where(relative_dist > 0.6, dye, 0)

    concentration       = np.zeros_like(local_radii)
    concentration_inner = np.zeros_like(local_radii)
    concentration_outer = np.zeros_like(local_radii)

    rs = local_radii[np.nonzero(local_radii)]
    coords_skel = np.transpose(np.nonzero(local_radii))

    for i in np.arange(0, len(rs), 1):
        y, x = coords_skel[i][0], coords_skel[i][1]
        disk = ((cx[np.newaxis,:]-x)**2 + (cy[:,np.newaxis]-y)**2 < rs[i]**2)

        concentration[y][x]       = np.sum(disk * dye)       / np.sum(disk)
        concentration_inner[y][x] = np.sum(disk * dye_inner) / np.sum(disk * dye_inner.astype(bool))
        concentration_outer[y][x] = np.sum(disk * dye_outer) / np.sum(disk * dye_outer.astype(bool))

    return concentration, concentration_inner, concentration_outer


def project_on_skeleton(dye, skeleton):
    ### orth. projection on skeleton by finding the nearest point in skeleton
    intensity   = np.zeros_like(dye)
    no          = np.zeros_like(dye)


    dist, inds = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                      return_distances=True, return_indices=True)

    dye_f = dye.flatten()
    inds0 = inds[0].flatten()
    inds1 = inds[1].flatten()

    for i in  range(0, len(dye_f)):
        if dye_f[i]!=0:
            intensity[inds0[i]][inds1[i]] += dye_f[i]
            no[inds0[i]][inds1[i]]+=1

    concentration = np.true_divide(intensity, no, out=np.zeros_like(no), where=intensity!=0)
    return concentration


def inter_mean(dye, skeleton, mask, local_radii, interval_size = 10):
    ### orth. projection on skeleton by finding the nearest point in skeleton
    # + average over nearest pixel in skeleton within r <= ~ 10 or so
    intensity_inner = np.zeros_like(dye)
    intensity_outer = np.zeros_like(dye)
    no_inner        = np.zeros_like(dye)
    no_outer        = np.zeros_like(dye)

    relative_dist = relative_distance(skeleton, mask, local_radii)

    inds = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                      return_distances=False, return_indices=True)

    dye_f = dye.flatten()
    inds0 = inds[0].flatten()
    inds1 = inds[1].flatten()
    relative_dist_f = relative_dist.flatten()

    for i in  range(0, len(dye_f)):
        if dye_f[i]!=0:
            if relative_dist_f[i]<0.6:
                intensity_inner[inds0[i]][inds1[i]] += dye_f[i]
                no_inner[inds0[i]][inds1[i]]+=1
            else:
                intensity_outer[inds0[i]][inds1[i]] += dye_f[i]
                no_outer[inds0[i]][inds1[i]]+=1


    intensity = intensity_inner + intensity_outer
    no = no_inner + no_outer

    concentration = np.true_divide(intensity, no,
                                   out=np.zeros_like(dye), where=intensity!=0)

    concentration_inner = np.true_divide(intensity_inner, no_inner,
                                         out=np.zeros_like(dye), where=intensity_inner!=0)

    concentration_outer = np.true_divide(intensity_outer, no_outer,
                                         out=np.zeros_like(dye), where=intensity_outer!=0)

    concentration         = disk_filter(concentration, concentration, interval_size) * skeleton
    concentration_inner   = disk_filter(concentration_inner, concentration, interval_size) * skeleton
    concentration_inner   = disk_filter(concentration_outer, concentration, interval_size) * skeleton


    return concentration, concentration_inner, concentration_outer
