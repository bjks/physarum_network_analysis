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

    kymos = [ [] for i in range(7) ]

    path                = []
    alignment           = []

    for set in data_sets:
        print(">> Analyse %s"  % (set.file_dat))

        values, seed_position, skeleton, l = extract_branch(set.file_dat,
                                                            seed_position,
                                                            ['local_radii',
                                                            'c_green',
                                                            'c_inner_green',
                                                            'c_outer_green',
                                                            'c_texas',
                                                            'c_inner_texas',
                                                            'c_outer_texas'])

        along_path , path_coords, last_endpoint = follow_all_paths(values, last_endpoint)
        path.append(path_coords)

        for i in range(7):
            kymos[i].append(along_path[i])

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
                        branch_map          = branch_map,
                        kymo_local_radii    = kymos[0],
                        kymo_c_green        = kymos[1],
                        kymo_inner_green    = kymos[2],
                        kymo_outer_green    = kymos[3],
                        kymo_c_texas        = kymos[4],
                        kymo_inner_texas    = kymos[5],
                        kymo_outer_texas    = kymos[6],
                        path                = path,
                        alignment           = alignment)


    ### plotting
    if not os.path.exists(set.path_plots):
        os.mkdir(set.path_plots)

    path_name = set.file_plot_set + '_branch_' + str(label) + '/'
    if not os.path.exists(path_name):
        os.mkdir(path_name)

    plot_branch(branch_map, label, path_name)


#########################################
############### skeleton ################
#########################################

def main(): ## python3 skeleton.py <seed index (0,...)>
    method = 'inter_mean'

    set_keyword = os.sys.argv[1]
    label       = int(os.sys.argv[2])
    
    seed_position = data(set_keyword).seed_positions[label]

    data_sets = [data(set_keyword, i, method, color='sep') for i in range(data(set_keyword).first, data(set_keyword).last)]
    process_skeleton(data_sets, seed_position, label)


if __name__ == "__main__":
    main()
