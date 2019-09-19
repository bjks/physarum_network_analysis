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
        if np.size(set.file_raw)>1:
            file_raw = file_raw[0]
        else:
            file_raw = set.file_raw
        network = invert_bf(read_file(file_raw))

    mask        = create_mask(network, set.sigma, set.threshold, set.halo_sig)
    mask        = extract_network(mask, set.extract)

    skeleton    = extract_skeleton(mask, method='medial_axis', branch_thresh=150)

    local_radii                    = extract_radii(mask, skeleton)

    rel_dist, radii_map            = relative_distance(skeleton, mask, local_radii)

    network_clean = np.multiply(network, mask)

    _, network_clean     = remove_spots(network_clean, mask, set.spots_radius, set.thresh_spots)


    if set.method == 'inter_mean':
        concentration, \
        concentration_inner, \
        concentration_outer = inter_mean(network_clean, skeleton, mask,
                                            local_radii, rel_dist,
                                            interval_size=50, div=0.5)



    map_names = np.array(['network_clean',
                        'skeleton',
                        'local_radii',
                        'mask',
                        'rel_dist',
                        'radii_map',
                        'concentration',
                        'concentration_inner',
                        'concentration_outer'])


    maps =      np.array([network_clean,
                        skeleton,
                        local_radii,
                        mask,
                        rel_dist,
                        radii_map,
                        concentration,
                        concentration_inner,
                        concentration_outer])


    if set.analyse_flow:
        bf_frames = [invert_bf(read_file(file)) for file in set.file_raw]
        flow_field_x, flow_field_y =  average_flow_over_frames(bf_frames, mask,
                                                upsample=1, window_size=30,
                                                sampling=20, search_extend=10)

        map_names = np.append(map_names, ['flow_field_x', 'flow_field_y'])
        maps      = np.append(maps, [flow_field_x, flow_field_y])


    to_save_dict = {n:m for n,m in zip(map_names, maps)}

    np.savez_compressed(set.file_dat,  **to_save_dict)



#########################################
##############   network   ##############
#########################################

def main(): ## python3 ratiometric.py <keyword> <first> <last(+1)>

    set_keyword = os.sys.argv[1]
    method = 'inter_mean'
    color =  str(os.sys.argv[2])

    data_sets = [data(set_keyword, i, method, color=color) for i in range( int(os.sys.argv[3]),int(os.sys.argv[4]) )]

    for set in data_sets:
        process_network(set)

if __name__ == "__main__":
    main()
