from analysis.network_analysis import *
from analysis.data_sets import *
import skimage.morphology as morph
import os
from analysis.plotting import *
from multiprocessing.dummy import Pool as ThreadPool


def process_ratiometric(set, use_ref_for_mask = True):
    ### read ###
    green = read_file(set.file_raw1)
    texas = read_file(set.file_raw2)

    ### mask, spots removal ###
    if use_ref_for_mask:
        mask = create_mask(texas, set.sigma, set.threshold, set.halo_sig)
    else:
        mask = create_mask(green, set.sigma, set.threshold, set.halo_sig)

    mask = extract_nerwork(mask)

    green_clean     = np.multiply(green, mask)
    texas_clean     = np.multiply(texas, mask)
    spots_mask, green_spotless     = remove_spots(green_clean, mask, set.spots_sig, set.thresh_spots)

    ### skeleton, radii ###
    skeleton                       = extract_skeleton(mask)
    local_radii                    = extract_radii(mask, skeleton)

    ratio                          = calc_ratio(green_spotless, texas_clean)

    ### projection methods ###
    if set.method == 'disk_mean':
        concentration, concentration_inner, concentration_outer = circle_mean(ratio, skeleton, mask, local_radii)

    if set.method == 'inter_mean':
        concentration, concentration_inner, concentration_outer = inter_mean(ratio, skeleton, mask, local_radii)

    np.savez_compressed(set.file_dat,           green_clean         = green_clean,
                                                texas_clean         = texas_clean,
                                                skeleton            = skeleton,
                                                local_radii         = local_radii,
                                                mask                = mask,
                                                green_spotless      = green_spotless,
                                                ratio               = ratio,
                                                concentration       = concentration,
                                                concentration_inner = concentration_inner,
                                                concentration_outer = concentration_outer)



#########################################
############## ratiometric ##############
#########################################

def main(): ## python3 ratiometric.py <keyword> <first> <last(+1)>

    set_keyword = os.sys.argv[1]
    method  = 'inter_mean'

##### green(i) texas(i+0.5), texas(i-0.5) green(i) ####
    for order in ['gt', 'tg']:
        data_sets = [data(set_keyword, i, method, color=order) for i in range( int(os.sys.argv[2]),int(os.sys.argv[3]) )]

        if not os.path.exists(data(set_keyword).path_results):
            os.mkdir(data(set_keyword).path_results)

        for set in data_sets:
            print(set.file_dat)
            process_ratiometric(set)


if __name__ == "__main__":
    main()
