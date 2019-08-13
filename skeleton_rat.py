from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *
import os

from multiprocessing.dummy import Pool as ThreadPool
import itertools


################## skeleton ####################
def process_skeleton(data_sets, seed_position, label):
    last_endpoint = seed_position

    kymograph_local_radii   = []
    kymograph_concentration = []
    kymograph_inner         = []
    kymograph_outer         = []

    path                    = []
    alignment               = []

    for set in data_sets:
        print(">> Analyse %s"  % (set.file_dat))

        values, seed_position, skeleton, l = extract_branch(set.file_dat, seed_position,
                                                            ['local_radii',
                                                            'concentration',
                                                            'concentration_inner',
                                                            'concentration_outer'])

        along_path , path_coords, last_endpoint = follow_all_paths(values, last_endpoint)
        path.append(path_coords)

        kymograph_local_radii.append(along_path[0])
        kymograph_concentration.append(along_path[1])
        kymograph_inner.append(along_path[2])
        kymograph_outer.append(along_path[3])
        ##### use first set as reference point: ####
        if set == data_sets[0]:
            alignment.append(int(len(path_coords)/2))
            reference_point = path_coords[int(len(path_coords)/2)]
            lenght0 = len(path_coords)
            branch_map = skeleton + values[0]

        else:
            point_to_align = closest_point_in_skel(reference_point, values[0])
            alignment.append(path_coords.index([point_to_align[0], point_to_align[1]]))
            if len(path_coords) < 0.2*lenght0:
                break

    ### save everything
    np.savez_compressed(set.file_dat_set + '_branch_' + str(label),
                        branch_map              = branch_map,
                        kymograph_local_radii   = kymograph_local_radii,
                        kymograph_concentration = kymograph_concentration,
                        kymograph_inner         = kymograph_inner,
                        kymograph_outer         = kymograph_outer,
                        path                    = path,
                        alignment               = alignment)


    ### plotting
    if not os.path.exists(set.path_plots):
        os.mkdir(set.path_plots)

    path_name = set.file_plot_set + '_branch_' + str(label) + '/'
    if not os.path.exists(path_name):
        os.mkdir(path_name)

    plot_branch(branch_map, label, path_name)

    plot_kymos(kymograph_local_radii, kymograph_concentration, 'radii', 'concentration',
               label, path_name, 'reference_point', alignment)

    plot_kymos(kymograph_inner, kymograph_outer, 'innner', 'outer',
               label, path_name, 'reference_point', alignment)


#########################################
############### skeleton ################
#########################################

def main(): ## python3 skeleton.py <keyword> <seed index (0,...)>
    method = 'inter_mean'

    set_keyword = os.sys.argv[1]
    color =   str(os.sys.argv[2])
    label =   int(os.sys.argv[3])
    seed_position = data(set_keyword).seed_positions[label]
    if color == 'both':
        colors = ['tg', 'gt']
    else:
        colors = [color]

    for color in colors:
        data_sets = [data(set_keyword, i, method, color=color) for i in range(data(set_keyword).first, data(set_keyword).last)]

        process_skeleton(data_sets, seed_position, label)


if __name__ == "__main__":
    main()
