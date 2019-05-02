from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.data_setsBF import *
import skimage.morphology as morph
import os
from analysis.plotting import *
from multiprocessing.dummy import Pool as ThreadPool

def invert_BF(image):
    return - image + np.max(image)


def process_network(set, color = 'texas'):

    if color=='texas':
        network = read_file(set.file_raw2)
    elif color=='green':
        network = read_file(set.file_raw1)
    elif color=='BF':
        network = invert_BF(read_file(set.file_raw))


    mask = create_mask(network, set.sigma, set.threshold, set.halo_sig )
    mask = extract_nerwork(mask)

    network_clean = np.multiply(network, mask)

    skeleton                       = extract_skeleton(mask)
    local_radii                    = extract_radii(mask, skeleton)


    if set.method == 'disk_mean':
        concentration, concentration_inner, concentration_outer = circle_mean(network_clean, skeleton, mask, local_radii)

    if set.method == 'inter_mean':
        concentration, concentration_inner, concentration_outer = inter_mean(network_clean, skeleton, mask, local_radii)

    np.savez_compressed(set.file_dat+color,     network_clean       = network_clean,
                                                skeleton            = skeleton,
                                                local_radii         = local_radii,
                                                concentration       = concentration,
                                                concentration_inner = concentration_inner,
                                                concentration_outer = concentration_outer)



#########################################
##############   network   ##############
#########################################

def main(): ## python3 ratiometric.py <keyword> <first> <last(+1)>

    set_keyword = os.sys.argv[1]
    method = 'inter_mean'
    color = 'BF'
    if color == 'BF':
        data_sets = [dataBF(set_keyword, i, method) for i in range( int(os.sys.argv[2]),int(os.sys.argv[3]) )]

        if not os.path.exists(dataBF(set_keyword).path_results):
            os.mkdir(dataBF(set_keyword).path_results)
    else:
        data_sets = [data(set_keyword, i, method) for i in range( int(os.sys.argv[2]),int(os.sys.argv[3]) )]

        if not os.path.exists(data(set_keyword).path_results):
            os.mkdir(data(set_keyword).path_results)

    for set in data_sets:
        process_network(set, color=color)

    print('Done')

if __name__ == "__main__":
    main()
