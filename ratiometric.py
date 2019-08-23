from analysis.network_analysis import *
from analysis.data_sets import *
import skimage.morphology as morph
import os
from analysis.plotting import *
from multiprocessing.dummy import Pool as ThreadPool


def process_ratiometric(set, use_ref_for_mask = True):
    print("ratiometric: ", set.file_dat)
    green = read_file(set.file_raw1)
    texas = read_file(set.file_raw2)

    green = background_correction(green, set.file_raw, set.sigma,
                                  set.lower_thresh, set.halo_sig)

    texas = background_correction(texas, set.file_raw, set.sigma,
                                  set.lower_thresh, set.halo_sig)

    if use_ref_for_mask:
        mask = create_mask(texas, set.sigma, set.threshold, set.halo_sig)
    else:
        mask = create_mask(green, set.sigma, set.threshold, set.halo_sig)

    # keep only the largest objects, number given by 'extract'
    mask = extract_network(mask, set.extract)

    ### skeleton, radii ### add branch_thresh to config files?!
    skeleton = extract_skeleton(mask, method='skeletonize', branch_thresh=500)


    local_radii = extract_radii(mask, skeleton)

    rel_dist, radii_map = relative_distance(skeleton, mask, local_radii)


    green_clean = np.multiply(green, mask)
    texas_clean = np.multiply(texas, mask)

    # _, green_clean = remove_spots(green_clean, mask, set.spots_radius, set.thresh_spots)
    # _, texas_clean = remove_spots(texas_clean, mask, set.spots_radius, set.thresh_spots)


    ratio                          = calc_ratio(green_clean, texas_clean)

    ### projection methods ###

    if set.method == 'inter_mean':
        c_green, \
        c_inner_green, \
        c_outer_green = inter_mean(green_clean, skeleton, mask, local_radii,
                                            rel_dist, interval_size=50, div=0.5)

        c_texas, \
        c_inner_texas, \
        c_outer_texas = inter_mean(texas_clean, skeleton, mask, local_radii,
                                            rel_dist, interval_size=50, div=0.5)


    np.savez_compressed(set.file_dat,   green_clean         = green_clean,
                                        texas_clean         = texas_clean,
                                        skeleton            = skeleton,
                                        local_radii         = local_radii,
                                        mask                = mask,
                                        ratio               = ratio,
                                        rel_dist            = rel_dist,
                                        radii_map           = radii_map,
                                        c_green             = c_green,
                                        c_inner_green       = c_inner_green,
                                        c_outer_green       = c_outer_green,
                                        c_texas             = c_texas,
                                        c_inner_texas       = c_inner_texas,
                                        c_outer_texas       = c_outer_texas)



#########################################
############## ratiometric ##############
#########################################

def main(): ## python3 ratiometric.py <keyword> <first> <last(+1)>

    set_keyword = os.sys.argv[1]
    start       = int(os.sys.argv[2])
    stop        = int(os.sys.argv[3])
    method  = 'inter_mean'

    data_sets = [data(set_keyword, i, method, color='sep') for i in range( start, stop )]

    for set in data_sets:
        process_ratiometric(set)

if __name__ == "__main__":
    main()
