import skimage.morphology as morph
import os
from multiprocessing.dummy import Pool as ThreadPool

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.flow_analysis import *
from analysis.plotting import *


def process_network(set):

    if set.color=='texas':
        network = read_file(set.file_raw2)
    elif set.color=='green':
        network = read_file(set.file_raw1)
    elif set.color=='bf':
        if np.size(set.file_raw)>1:
            file_raw = set.file_raw[0]
        else:
            file_raw = set.file_raw
        network = invert_bf(read_file(file_raw))


    ### keep only the n largest objects, n given by 'extract'
    mask = create_mask(network, set.sigma, set.threshold, set.halo_sig)
    mask = extract_network(mask, set.extract)

    ### skeleton, radii ###
    skeleton = extract_skeleton(mask, method='skeletonize',
                                branch_thresh = set.branch_thresh,
                                extract = set.extract)

    local_radii                    = extract_radii(mask, skeleton)
    rel_dist, radii_map            = relative_distance(skeleton, mask, local_radii)

    network_clean = np.multiply(network, mask)

    _, network_clean     = remove_spots(network_clean, mask, set.spots_radius, set.thresh_spots)


    if set.method == 'inter_mean':
        concentration, \
        concentration_inner, \
        concentration_outer = inter_mean(network_clean, skeleton, mask, rel_dist,
                                            interval_size=50, div=0.5,
                                            corr_for_missing_branches = False)


        # gathering data to create dictionary for np.savez, data is not copied!
        maps = [('set', set),
                    ('network_clean', network_clean),
                    ('skeleton', skeleton),
                    ('local_radii', local_radii),
                    ('mask', mask),
                    ('rel_dist', rel_dist),
                    ('radii_map', radii_map),
                    ('concentration', concentration),
                    ('concentration_inner', concentration_inner),
                    ('concentration_outer', concentration_outer)]


    if set.analyse_flow:
        bf_frames = [invert_bf(read_file(file)) for file in set.file_raw]
        flow_field_x, flow_field_y = average_flow_over_frames(bf_frames, mask,
                                                            upsample=1,
                                                            window_size=30,
                                                            sampling=10,
                                                            search_extend=20,
                                                            corr_tresh=0.6,
                                                            outlier_thresh=None)

        flow_x = inter_mean(flow_field_x, skeleton, mask, rel_dist,
                                interval_size=50, div=0.5,
                                corr_for_missing_branches = True,
                                return_only_conc = True)
        flow_y = inter_mean(flow_field_y, skeleton, mask, rel_dist,
                                interval_size=50, div=0.5,
                                corr_for_missing_branches = True,
                                return_only_conc = True)


        maps.extend((('flow_x', flow_x), ('flow_y', flow_y)))


    ######## SAVE #######
    to_save_dict = dict(maps)
    np.savez_compressed(set.file_dat,  **to_save_dict)



#########################################
##############   network   ##############
#########################################

def main(): ## python3 ratiometric.py <keyword> <first> <last(+1)>

    set_keyword = os.sys.argv[1].strip()
    color       = os.sys.argv[2].strip()
    start       = int(os.sys.argv[3].strip())
    stop        = int(os.sys.argv[4].strip())
    method = 'inter_mean'


    data_sets = [data(set_keyword, i, method, color=color) for i in range(start, stop)]

    for set in data_sets:
        process_network(set)



if __name__ == "__main__":
    main()
