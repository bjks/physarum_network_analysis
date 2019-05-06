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
    set_keyword = os.sys.argv[1]
    method = 'inter_mean'
    label = int(os.sys.argv[2])

    orders = ['gt', 'tg']
    seed_position = data(set_keyword).seed_positions[label]
    for order in orders:
        data_sets = [data(set_keyword, i, method, color=order) for i in range(data(set_keyword).first, data(set_keyword).last)]
        process_skeleton(data_sets, seed_position, label)

if __name__ == "__main__":
    main()
