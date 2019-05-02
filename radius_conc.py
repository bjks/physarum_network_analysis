from analysis.data_sets import *
from analysis.tools import *
from analysis.network_analysis import *
from analysis.skeleton_analysis import *
from analysis.plotting import *
import os

from multiprocessing.dummy import Pool as ThreadPool
from scipy.optimize import curve_fit
import itertools

def lin_regr(x, m, c):
    return m*x + c


################## skeleton ####################
def collect_data(data_sets, color, keyword, seed_position):

    data = np.array([])

    for set in data_sets:
        print("\r>> Analyse %s"  % (set.file_dat), end='')
        set_data   = np.load(set.file_dat + color + '.npz')[keyword]
        skeleton   = np.load(set.file_dat + color + '.npz')['skeleton']
        cx = np.arange(0, skeleton.shape[1])
        cy = np.arange(0, skeleton.shape[0])

        r = 500
        y = 1300
        x = 900

        roi = ((cx[np.newaxis,:]-x)**2 + (cy[:,np.newaxis]-y)**2 < r**2)

        seed_position   = closest_point_in_skel(seed_position, skeleton)
        branches        = label_branches(skeleton)
        label           = branches[seed_position[0], seed_position[1]]
        selected_branch = np.where((branches==label) , 1, 0) * roi

        list = set_data[np.nonzero(set_data*selected_branch)]
        data = np.append(data, list)

    return data



#########################################
############### skeleton ################
#########################################

def main(): ## python3 skeleton.py <keyword> <seed index (0,...)>
    set_keyword = os.sys.argv[1]
    method = 'projection'

    data_sets = [data(set_keyword, i, method) for i in range(data(set_keyword).first, 10)]
    # data_sets = [data(set_keyword, i, method) for i in range(1,2)]
    seed_position = data(set_keyword).seed_positions[int(os.sys.argv[2])]

    c = collect_data(data_sets, 'texas', 'concentration', seed_position)
    r = collect_data(data_sets, 'texas', 'local_radii', seed_position)

    plt.scatter(c, r, s=0.3)

    scaling, offset = curve_fit(lin_regr, c, r)[0]

    print(scaling)
    plt.show()

    print('Done')


if __name__ == "__main__":
    main()
