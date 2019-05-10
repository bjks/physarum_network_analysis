from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *
import os

from multiprocessing.dummy import Pool as ThreadPool
import itertools


#########################################
############### skeleton ################
#########################################

def main(): ## python3 skeleton.py <keyword> <seed index (0,...)>
    method = 'inter_mean'

    set_keyword = os.sys.argv[1]
    color =  str(os.sys.argv[2])
    label = int(os.sys.argv[3])
    seed_position = data(set_keyword).seed_positions[label]

    if color == 'both':
        for order in ['gt', 'tg']:
            data_sets = [data(set_keyword, i, method, color=order) for i in range(data(set_keyword).first, data(set_keyword).last)]
            process_skeleton(data_sets, seed_position, label)
    else:
        data_sets = [data(set_keyword, i, method, color=color) for i in range(data(set_keyword).first, data(set_keyword).last)]
        process_skeleton(data_sets, seed_position, label)



if __name__ == "__main__":
    main()
