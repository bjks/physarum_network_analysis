import skimage.morphology as morph
import os

import sys
sys.path.append("..")

from analysis.network_analysis import *
from analysis.data_sets import *
from analysis.flow_analysis import *
from analysis.tools import *

def fraction_unmasked(arr, fraction=0.1):
    avg = np.mean(arr)
    if avg > fraction:
        return True
    else:
        return False

def process_network(set, set2):

    if not os.path.isfile(set.file_dat + '.npz'):
        # 'frame' is the refernce frame, i.e if 'frame' corresponds to t=1
        # everything is saved in *t0001*
        frame = invert_bf(read_file(set.file_raw))
        frame2 = invert_bf(read_file(set2.file_raw))

        ### keep only the n largest objects, n given by 'extract'
        mask = create_mask(frame, set.sigma, set.threshold, set.halo_sig)
        mask = extract_network(mask, set.extract)

        flow_x, flow_y = flow_quantification(frame, frame2, mask,
                                                        upsample=1,
                                                        window_size=30,
                                                        sampling=10,
                                                        search_extend=20,
                                                        return_scalar_v=False,
                                                        corr_tresh=0.3,
                                                        outlier_thresh=None,
                                                        mask_crit=np.any)



        maps = [('mask', mask.astype(bool)),
                ('flow_x', flow_x),
                ('flow_y', flow_y) ]


        ######## SAVE #######
        to_save_dict = dict(maps)
        np.savez_compressed(set.file_dat,  **to_save_dict)



#########################################
##############   network   ##############
#########################################

def main(): ## python3 piv.py <keyword> <color> <first> <last(+1)>

    print(os.sys.argv)
    set_keyword = os.sys.argv[1].strip()
    color       = os.sys.argv[2].strip()
    start       = int(os.sys.argv[3].strip())
    stop        = int(os.sys.argv[4].strip())

    for i in range(start, stop):
        process_network(data(set_keyword, i,    method='piv', color=color),
                        data(set_keyword, i+1,  method='piv', color=color) )

    log_message(data_sets[0], 'piv', start, stop)


if __name__ == "__main__":
    main()
