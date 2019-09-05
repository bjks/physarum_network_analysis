from analysis.network_analysis import *
from analysis.data_sets import *
import skimage.morphology as morph
import os
from analysis.plotting import *
from multiprocessing.dummy import Pool as ThreadPool


def process_network(set):

    if set.color=='texas':
        network = read_file(set.file_raw2)
    elif set.color=='green':
        network = read_file(set.file_raw1)
    elif set.color=='bf':
        network = invert_bf(read_file(set.file_raw))

    mask        = create_mask(network, set.sigma, set.threshold, set.halo_sig)
    mask        = extract_network(mask, set.extract)

    skeleton    = extract_skeleton(mask, method='medial_axis', branch_thresh=150)

    local_radii                    = extract_radii(mask, skeleton)

    rel_dist, radii_map            = relative_distance(skeleton, mask, local_radii)

    network_clean = np.multiply(network, mask)

    _, network_clean     = remove_spots(network_clean, mask, set.spots_radius, set.thresh_spots)



    if set.method == 'disk_mean':
        concentration, \
        concentration_inner, \
        concentration_outer = circle_mean(network_clean, skeleton, mask,
                                            local_radii, rel_dist, div=0.5)

    if set.method == 'inter_mean':
        concentration, \
        concentration_inner, \
        concentration_outer = inter_mean(network_clean, skeleton, mask,
                                            local_radii, rel_dist,
                                            interval_size=50, div=0.5)

    np.savez_compressed(set.file_dat,   network_clean       = network_clean,
                                        skeleton            = skeleton,
                                        local_radii         = local_radii,
                                        mask                = mask,
                                        rel_dist            = rel_dist,
                                        radii_map           = radii_map,
                                        concentration       = concentration,
                                        concentration_inner = concentration_inner,
                                        concentration_outer = concentration_outer)



#########################################
##############   network   ##############
#########################################

def main(): ## python3 ratiometric.py <keyword> <first> <last(+1)>

    set_keyword = os.sys.argv[1]
    method = 'inter_mean'
    color =  str(os.sys.argv[2])

    data_sets = [data(set_keyword, i, method, color=color) for i in range( int(os.sys.argv[3]),int(os.sys.argv[4]) )]

    for set in data_sets:
        print(set.file_dat)
        process_network(set)

    print('Done')

if __name__ == "__main__":
    main()
