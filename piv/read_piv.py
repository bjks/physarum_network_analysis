import skimage.morphology as morph
import os

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.flow_analysis import *

def analyse_flow_field(data_sets):

    fields_x = np.array([np.load(set.file_dat + '.npz')['flow_x'] for set in data_sets])
    fields_y = np.array([np.load(set.file_dat + '.npz')['flow_y'] for set in data_sets])
    masks =  np.array([np.load(set.file_dat + '.npz')['mask'] for set in data_sets])

    show_im(masks[0])

    mean_x, mean_y = [], []
    for i in range(int(len(fields_x)/10)):
        mean_x.append(np.sum(fields_x[i:i+10], axis=0))
        mean_y.append(np.sum(fields_y[i:i+10], axis=0))


    for my, mx in zip(mean_y, mean_x):
        abs = (my + mx)**2
        plt.imshow( abs , vmax=np.mean(abs)+np.std(abs), vmin=np.mean(abs)-np.std(abs) )
        plt.show()


#########################################
##############   network   ##############
#########################################

def main(): ## python3 piv.py <keyword> <color> <first> <last(+1)>

    set_keyword = os.sys.argv[1].strip()
    color       = os.sys.argv[2].strip()

    indx = range(data(set_keyword).first, data(set_keyword).last)
    data_sets = [data(set_keyword, i, method='piv', color=color) for i in indx]

    analyse_flow_field(data_sets)


if __name__ == "__main__":
    main()
