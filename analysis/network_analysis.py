import numpy as np
import matplotlib.pyplot as plt

import skimage.morphology as morph              # skeleton
from skimage.measure import regionprops         # extract_network
from skimage.filters import threshold_otsu      # mask
from skimage.feature import peak_local_max


from scipy import ndimage as ndi                    # gaussian_filter
from astropy.convolution import convolve_fft

############ File reading ########
def greyscale(image):
    return np.dot(image[...,:3], [0.299, 0.587, 0.114])

def read_file(filename):
    image = plt.imread(filename, 'tif').astype(float)
    if image.ndim==3:
        return greyscale(image)
    else:
        return image

def invert_bf(image):
    return - image + np.max(image)

###### visualization tools ######
### create diluted skeleton and masks the background for plotting ###
def thick_skeleton(skeleton, times = 10):

    if times > 0:
        selem = morph.disk(times)
        thick_skeleton = morph.dilation(skeleton.astype(float), selem)

    else:
        thick_skeleton = skeleton.copy().astype(float)

    thick_skeleton = np.ma.masked_where(thick_skeleton == 0, thick_skeleton)
    return thick_skeleton
    

def show_im(image, skel=False):
    plt.close()
    if skel:
        plt.imshow(thick_skeleton(image))
    else:
        plt.imshow(image)
    plt.colorbar()
    plt.show()
    plt.close()

#####################################
############# Creating mask #########
#####################################
# mask consists of doubles to avoid stupid buck in the future,
# although it costs some storage...

def create_mask(dye, sig, thresh=None, halo_sig=None):
    ### if no threshold is given, thresh is determined using Otsu’s method
    if thresh == None:
        smooth = ndi.gaussian_filter(dye, sigma=sig)
        thresh = threshold_otsu(smooth)
        mask = np.where(smooth > thresh, 1., 0.)


    ### if halo_sig is given the local background is estimated and used for mask
    elif halo_sig != None:
        halo = ndi.gaussian_filter(dye, sigma=halo_sig)
        smooth = ndi.gaussian_filter(dye, sigma=sig)
        mask = np.where(smooth > thresh * halo, 1., 0.)


    ### pixels are compared to mean of entire image
    else:
        smooth = ndi.gaussian_filter(dye, sigma=sig)
        mask = np.where(smooth > thresh*(np.mean(smooth)), 1., 0.)

    return mask

def extract_network(mask, n):
    if n == None:
        return mask

    # returns mask that contains only the n largest connected areas
    labels = morph.label(mask, connectivity=2)
    regions = regionprops(labels)

    areas = np.argsort([r.area for r in regions])
    #list of labels of the n largest areas
    selected_labels = [regions[a].label for a in areas[-n:]]

    return np.where(np.isin(labels, selected_labels), 1., 0.)

#####################################
######## image interpolation ########
########### discontinued!! ##########
#####################################
def interpl_masks(mask1, mask2):
    # finds outline in the center between the outlines of two maskss
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

def node_detection(bitmap):
    skeleton = np.where(bitmap!=0, 1, 0)
    sum =   (np.roll(skeleton,  1, axis=1) +
             np.roll(skeleton, -1, axis=1) +
             np.roll(skeleton,  1, axis=0) +
             np.roll(skeleton, -1, axis=0) +
             np.roll(np.roll(skeleton, 1, axis=0), 1, axis=1) +
             np.roll(np.roll(skeleton, 1, axis=0), -1, axis=1) +
             np.roll(np.roll(skeleton, -1, axis=0), 1, axis=1) +
             np.roll(np.roll(skeleton, -1, axis=0), -1, axis=1))

    nodes       = np.where(sum > 2, 1, 0) * skeleton
    endpoints   = np.where(sum ==1, 1, 0) * skeleton

    return nodes, endpoints


def remove_branch(branch, endpoints, branch_thresh):
    if np.max(branch + endpoints) > 1 and np.sum(branch) < branch_thresh:
        return True
    else:
        return False

def set_borders_to(arr, value=0):
    new_array = arr.copy()
    new_array[0,:]=value
    new_array[-1,:]=value
    new_array[:,0]=value
    new_array[:,-1]=value
    return new_array

def reflect_boundaries(arr):
    lr = np.fliplr(arr)
    new = np.concatenate((lr, arr, lr), axis=1)
    ud = np.flipud(new)
    return np.concatenate((ud, new, ud), axis=0)

def inverse_reflect_boundaries(arr):
    return np.hsplit(np.vsplit(arr,3)[1], 3)[1]


def extract_skeleton(mask, method='medial_axis', branch_thresh=50,
                        extract=1, reflect=True):

    refl_mask = reflect_boundaries(mask)

    if method=='medial_axis':
        # returns the line consiting of pixels that have 2 (or more)
        # nearest pixels; often many small branches emerge
        medial_axis = morph.medial_axis(refl_mask)

    elif method=='skeletonize':
        # returns the skeleton calculated via morph. thinning, which does not
        # guarantee to get the center line in a pixel(!) image
        medial_axis = morph.skeletonize(refl_mask)

    medial_axis = inverse_reflect_boundaries(medial_axis)
    medial_axis = set_borders_to(medial_axis, value=0)

    if branch_thresh > 0:
        nodes, endpoints = node_detection(medial_axis)
        seperated_skel = medial_axis - nodes

        seper_l, no_l = morph.label(seperated_skel, connectivity=2, return_num=True)

        # iterate over branches
        for l in range(1, no_l+1):
            branch = np.where(seper_l==l, 1, 0)
            # remove branches that don't connect anything and are shorter than
            # branch_thresh
            if remove_branch(branch, endpoints, branch_thresh):
                medial_axis = np.where(branch==1, 0, medial_axis)

        # add nodes to conncect network again
        medial_axis = np.logical_or(medial_axis, nodes)
        # remove nodes that are not longer part of the network
        medial_axis = morph.remove_small_objects(medial_axis, 6, connectivity=2)

        # thinning necessary since nodes that lost a branch are wider that 1px
        medial_axis = morph.skeletonize(medial_axis)

    medial_axis = extract_network(medial_axis, extract)
    return medial_axis


def extract_radii(mask, skeleton):
    distance = ndi.distance_transform_edt(mask)
    local_radii = distance * skeleton # only keeps pixel in skeleton
    return local_radii

#########################################################################
########## replace spots with estim. background, based on nearby pixels #

def disk_filter(dye, mask, r):
    dye_filtered = np.where(mask!=0, dye, np.nan) # 0 -> NaN (to be ignored)

    kernel = morph.disk(r)
    dye_filtered = convolve_fft(dye_filtered, kernel, fill_value=np.nan,
                                quiet=True)

    return np.nan_to_num(dye_filtered) # NaN -> 0


def remove_spots(dye, mask, spots_sig, thresh_spots):
    #### calculate background ####
    background   = disk_filter(dye, mask, spots_sig)

    #### calculate mask withou spots ####
    spots_mask   = np.where(dye < thresh_spots * background, 1., 0.) * mask

    #### replace spots with estimated background ###
    estim_back   = disk_filter(dye, spots_mask, spots_sig)
    dye_spotless = np.where(spots_mask!=0, dye, estim_back ) * mask

    return spots_mask, dye_spotless

################################################
############ background_correction #############
################################################

def estimate_background(dye, mask, halo_sig):
    mask_inv    = np.where(mask==1., 0, 1)

    estim_back  = disk_filter(dye, mask_inv, halo_sig)
    return estim_back


def background_correction(dye, file_raw, sigma, lower_thresh, halo_sig):


    if np.size(file_raw)>1:
        file_raw = file_raw[0]

    bf          = read_file(file_raw)

    mask_bf     = create_mask(invert_bf(bf), sigma, lower_thresh, halo_sig)

    # estimate background of bf and dye
    back_bf     = estimate_background(bf, mask_bf, halo_sig)
    back_dye    = estimate_background(dye, mask_bf, halo_sig)

    # based on the assumption:
    # background (ligth below tube) / background that reaches the cam is constant
    # ie not dependend on the intensity and the wavelength
    # -> I_b,bf / I^0_b,bf * I^0_b, green = I_b,green
    added_back  = calc_ratio(bf, back_bf) * back_dye


    # tube signal = signal - background contribution
    corrected_dye = dye - added_back
    # necessary, otherwise negative values mess with mask
    # (not problematic, since it only concerns pixel outside network)
    return np.where(corrected_dye < 0, 0, corrected_dye)

#####################################
############ ratio calc #############
#####################################

def calc_ratio(dye, ref):
    # ratio between two arrays, returns zeros where denominator == 0
    return np.true_divide(dye, ref, out=np.zeros_like(dye), where=ref!=0)

#####################################
########## skeleton mapping #########
#####################################

def tube_radius_at_point(mask, skel, local_radii):
    ### calculate the radius (or other given quantity) at every position
    #   in network based on value of the nearest skeleton point
    radii = mask.copy().astype(float)
    inds = ndi.distance_transform_edt(np.invert(skel.astype(bool)),
                                      return_distances=False,
                                      return_indices=True)

    radii = radii.flatten()
    inds0 = inds[0].flatten()
    inds1 = inds[1].flatten()

    for i in  range(0, len(radii)):
        if radii[i]!=0:
            radii[i] = local_radii[inds0[i]][inds1[i]]

    return radii.reshape(np.shape(mask))

def relative_distance(skeleton, mask, local_radii):
    distance = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                          return_distances=True,
                                          return_indices=False)
    distance *= mask
    radii = tube_radius_at_point(mask, skeleton, local_radii)
    return np.true_divide(distance, radii, out=np.zeros_like(radii),
                            where=radii!=0), radii

def absolute_distance(skeleton, mask):
    distance = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                          return_distances=True,
                                          return_indices=False)
    return distance



def circle_mean(dye, skeleton, mask, local_radii, relative_dist, div=0.5):
    ###### 'sliding disk' #####
    cx = np.arange(0, dye.shape[1])
    cy = np.arange(0, dye.shape[0])

    #### calc realtive distance (r/R) to skeleton at every position
    # relative_dist, radii = relative_distance(skeleton, mask, local_radii)

    ### seperate outer and inner of network based on relative_dist
    dye_inner = np.where(relative_dist <= div, dye, 0)
    dye_outer = np.where(relative_dist > div, dye, 0)

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

    return  concentration, \
            concentration_inner, \
            concentration_outer



def project_on_skeleton(dye, skeleton):
    """ orth. projection on skeleton by finding the nearest point in skeleton
    only a few points are averaged, thus the results are noisy
    """

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

    concentration = np.true_divide(intensity, no, out=np.zeros_like(no),
                                   where=intensity!=0)
    return concentration



def inter_mean(dye, skeleton, mask, relative_dist,
                interval_size = 10, div = 0.5,
                corr_for_missing_branches = False,
                return_only_conc = False):
    """ orth. projection on skeleton by finding the nearest point in skeleton
     plus average over nearest pixel in skeleton within interval <= 10 pixel
     allow enough pixel to be considered such that "inner" "outer" can be
     returned.
     corr_for_missing_branches -> removes pixels that belong to branches that
     were removes
     """

    if corr_for_missing_branches:
         dye = np.where(relative_dist > 1, 0, dye)

    intensity_inner = np.zeros_like(dye)
    intensity_outer = np.zeros_like(dye)
    no_inner        = np.zeros_like(dye)
    no_outer        = np.zeros_like(dye)

    # relative_dist, radii = relative_distance(skeleton, mask, local_radii)

    inds = ndi.distance_transform_edt(np.invert(skeleton.astype(bool)),
                                      return_distances=False, return_indices=True)

    dye_f = dye.flatten()
    inds0 = inds[0].flatten()
    inds1 = inds[1].flatten()
    relative_dist_f = relative_dist.flatten()

    for i in  range(0, len(dye_f)):
        if dye_f[i]!=0:
            if relative_dist_f[i] < div:
                intensity_inner[inds0[i]][inds1[i]] += dye_f[i]
                no_inner[inds0[i]][inds1[i]]+=1
            else:
                intensity_outer[inds0[i]][inds1[i]] += dye_f[i]
                no_outer[inds0[i]][inds1[i]]+=1


    intensity = intensity_inner + intensity_outer
    no = no_inner + no_outer

    concentration = np.true_divide(intensity, no,
                                   out=np.zeros_like(dye),
                                   where=intensity!=0)

    concentration_inner = np.true_divide(intensity_inner, no_inner,
                                         out=np.zeros_like(dye),
                                         where=intensity_inner!=0)

    concentration_outer = np.true_divide(intensity_outer, no_outer,
                                         out=np.zeros_like(dye),
                                         where=intensity_outer!=0)

    if return_only_conc:
        return disk_filter(concentration, skeleton, interval_size) * skeleton


    concentration         = disk_filter(concentration,
                                        skeleton, interval_size) * skeleton
    concentration_inner   = disk_filter(concentration_inner,
                                        skeleton, interval_size) * skeleton
    concentration_outer   = disk_filter(concentration_outer,
                                        skeleton, interval_size) * skeleton


    return  concentration, \
            concentration_inner, \
            concentration_outer
