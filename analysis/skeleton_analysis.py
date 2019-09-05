import numpy as np
import skimage.morphology as morph
import copy
from analysis.network_analysis import *

from analysis.plotting import *
import os

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
    # finds endpoint of branch that is closest to the point 'last_endpoint'
    _, ends = node_detection(branch)
    coords_ends = np.transpose(np.nonzero(ends))
    dists = np.sum((coords_ends-last_endpoint)**2, axis=1)
    return coords_ends[np.argmin(dists)]


def follow_path(selected_branch, last_endpoint):
    # relurns values along the branch
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
    # relurns values along the path given by selected_branches
    # for all values in eg. selected_branches = [color, pressure, temperatur]
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

    #### replace selected_branch with shortest_path between two selected endpoints?
    # match 2 endpoints instead of seed_position ? (seed_position -> endpoints (2,2) )

    values_in_branch = []
    for q in quantities:
        values_in_branch.append(np.load(file_dat + '.npz')[q].astype(float) * selected_branch)

    return values_in_branch, seed_position, skeleton, label
