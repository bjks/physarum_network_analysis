import numpy as np
import skimage.morphology as morph
import copy
from analysis.plotting import *
import os


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


def closest_point_in_skel(seed, skeleton):
    coords_skel = np.transpose(np.nonzero(skeleton))
    dists = np.sum((coords_skel - seed)**2, axis=1)
    return coords_skel[np.argmin(dists)]


def label_branches(skeleton):
    nodes,_ = node_detection(skeleton)
    seperated_branches = skeleton - nodes
    return morph.label(seperated_branches, connectivity=2)


def next_pixel(seed, skel):
    i, j = seed[0], seed[1]
    directions =  [[i+1,j+1],
                    [i+1,j-1],
                    [i+1,j],
                    [i-1,j+1],
                    [i-1,j-1],
                    [i-1,j],
                    [i,j+1],
                    [i,j-1]     ]
    neigbours = ([skel[i[0],i[1]] for i in directions])
    x = np.nonzero(neigbours)

    if np.shape(x)[1]!=0:
        return directions[x[0][0]]
    else:
        return None

def find_branch_endpoints(branch, last_endpoint):
    _, ends = node_detection(branch)
    coords_ends = np.transpose(np.nonzero(ends))
    dists = np.sum((coords_ends-last_endpoint)**2, axis=1)
    return coords_ends[np.argmin(dists)]


def follow_path(selected_branch, last_endpoint):
    branch = selected_branch.copy()
    path_coords = []
    path_vals   = []
    endpoint = find_branch_endpoints(branch, last_endpoint)
    seed = last_endpoint = [endpoint[0], endpoint[1]]

    while seed!=None:
        path_coords.append(seed)
        path_vals.append(branch[seed[0],seed[1]])
        branch[seed[0],seed[1]] = 0
        seed = next_pixel(seed, branch)
    return path_vals, path_coords, last_endpoint

def follow_all_paths(selected_branches, last_endpoint):
    branches = copy.deepcopy(selected_branches)
    path_coords = []
    path_vals   = [ [] for b in branches ]
    endpoint = find_branch_endpoints(branches[0], last_endpoint)
    seed = last_endpoint = [endpoint[0], endpoint[1]]

    while seed!=None:
        path_coords.append(seed)

        for values, branch in zip(path_vals, branches):
            values.append(branch[seed[0],seed[1]])

        branches[0][seed[0],seed[1]] = 0
        seed = next_pixel(seed, branches[0])
    return path_vals, path_coords, last_endpoint


def extract_branch(file_dat, seed_position, quantities):
    # selects branch based on seed_position
    skeleton        = np.load(file_dat + '.npz')['skeleton'].astype(int)
    seed_position   = closest_point_in_skel(seed_position, skeleton)
    branches        = label_branches(skeleton)
    label           = branches[seed_position[0], seed_position[1]]
    selected_branch = np.where(branches==label, 1, 0)

    values_in_branch = []
    for q in quantities:
        values_in_branch.append(np.load(file_dat + '.npz')[q].astype(float) * selected_branch)

    return values_in_branch, seed_position, skeleton, label


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
