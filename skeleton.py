from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *
import os

import multiprocessing
import itertools


def quantities_to_analyse(set, prefix=''):

    quantities = np.array(['local_radii'])

    if set.color=='bf' or set.color=='texas' or set.color=='green':
        quantities = np.append(quantities,['concentration',
                                            'concentration_inner',
                                            'concentration_outer'])

    else:
        quantities = np.append(quantities,[ 'c_green',
                                            'c_inner_green',
                                            'c_outer_green',
                                            'c_texas',
                                            'c_inner_texas',
                                            'c_outer_texas'])

    if set.analyse_flow:
        quantities = np.append(quantities, ['flow_x', 'flow_y'])

    quantities = [prefix + q for q in quantities]

    return quantities, range(len(quantities))




################## skeleton ####################
def process_skeleton(data_sets, seed_position, label):
    last_endpoint = seed_position

    quantities, q_range = quantities_to_analyse(data_sets[0])

    kymo_names, _ =  quantities_to_analyse(data_sets[0], prefix='kymo_')

    kymos = [ [] for i in quantities]
    path                = []
    alignment           = []

    for set in data_sets:
        print(">> Analyse" + set.file_dat, label )

        values, seed_position, skeleton, l = extract_branch(set.file_dat,
                                                            seed_position,
                                                            quantities)

        along_path , path_coords, last_endpoint = follow_all_paths(values, last_endpoint)
        path.append(path_coords)

        for i in q_range:
            kymos[i].append(along_path[i])

        ##### use first set as reference point: ####
        if set == data_sets[0]:
            alignment.append(int(len(path_coords)/2))
            reference_point = path_coords[int(len(path_coords)/2)]
            lenght0 = len(path_coords)
            branch_map = skeleton.astype(int) + np.where(values[0]>0, 1, 0)

        else:
            point_to_align = closest_point_in_skel(reference_point, values[0])
            alignment.append(path_coords.index([point_to_align[0], point_to_align[1]]))
            if len(path_coords) < 0.2*lenght0:
                break

    ### save everything
    to_save_dict = {n:k for n,k in zip(kymo_names, kymos)}

    np.savez_compressed(branch_datfile(set, label),
                        branch_map          = branch_map,
                        path                = path,
                        alignment           = alignment,
                        **to_save_dict)


    ### plotting
    if not os.path.exists(set.path_plots):
        os.mkdir(set.path_plots)

    path_name = branch_plotpath(set, label)


    plot_branch(branch_map, label, path_name)


#########################################
############### skeleton ################
#########################################

def main(): ## python3 skeleton.py NAME
    method = 'inter_mean'

    set_keyword = os.sys.argv[1].strip()
    # label       = int(os.sys.argv[2])
    color       = 'sep'

    data_sets = [data(set_keyword, i, method, color=color) for i in range(data(set_keyword).first, data(set_keyword).last)]
    seed_positions, labels = get_seeds_positions(data_sets[0])

    # seed_position = data(set_keyword).seed_positions[label]
    # process_skeleton(data_sets, seed_positions, label)

    num_threads = multiprocessing.cpu_count()
    print("Number of detected cores: ", num_threads)

    p = multiprocessing.Pool(num_threads)
    p.starmap(process_skeleton, zip(itertools.repeat(data_sets),
                                    seed_positions,
                                    labels))
    p.close()
    p.join()


if __name__ == "__main__":
    main()
