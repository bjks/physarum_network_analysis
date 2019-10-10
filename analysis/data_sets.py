import sys
import os
import numpy as np


############### Tools ###################
def mk_mising_dir(path_name):
    if not os.path.exists(path_name):
        os.mkdir(path_name)
    return path_name

#####################################################
################ Data Structure #####################
#####################################################



class data_paths:
    def __init__(self, image_prefix, no, zeros, texas, green, bf):

        self.core1           = image_prefix + str(no).zfill(zeros)
        if self.color == 'gt':
            self.core2           = image_prefix + str(no+1).zfill(zeros)
        else:
            self.core2           = image_prefix + str(no).zfill(zeros)

        if self.subdir == None:
            dir = self.path
        else:
            dir = self.path + '/' + self.subdir

        if sys.platform.startswith('linux'): # i.e. on the (ubuntu) cluster
            dir_raw             = '/data.bpm/zeiss/BjoernZ/' + dir + '/'
            dir_results         = '/data.bpm/bksche/results/'
            dir_plots           = '/data.bpm/bksche/plots/'

        else:                               # local
            dir_raw             = '../image_data/'     + dir + '/'
            dir_results         = '../results/'
            dir_plots           = '../plots/'


        for d in [dir_results, dir_plots]:
            mk_mising_dir(d + self.path)
            if self.subdir != None:
                mk_mising_dir(d + self.path + '/' +self.subdir)


        if bf != None:
            bf = [x.strip() for x in bf.split()]
            if np.size(bf)>1:
                self.file_raw = [dir_raw + self.core1 + b + '.tif' for b in bf]
            elif np.size(bf)==1:
                self.file_raw = dir_raw + self.core1 + bf[0] + '.tif'
            else:
                self.file_raw = dir_raw + self.core1 + '.tif'

        if green != None:
            self.file_raw1      = dir_raw + self.core1 + green + '.tif'
        if texas != None:
            self.file_raw2      = dir_raw + self.core2 + texas + '.tif'





        self.path_results   = '../results/'    + dir + '/'
        self.path_plots     = '../plots/'      + dir + '/'

        full_file = self.core1   + '_' + self.method + '_' + self.color
        set_file  = 'set'        + '_' + self.method + '_' + self.color

        if self.lower_thresh != None:
            full_file += '_back_corr'
            set_file  += '_back_corr'

        self.file_dat       = self.path_results + full_file
        self.file_dat_flow  = self.path_results + full_file + '_flow'

        self.file_dat_set   = self.path_results + set_file

        self.file_plot              = mk_mising_dir(self.path_plots + 'samples/') + full_file
        self.file_plot_tube_profile = mk_mising_dir(self.path_plots + 'tube_profile/') + set_file
        self.file_plot_set          = self.path_plots + set_file



#####################################################
###################### data sets ####################
#####################################################
import configparser
import json

class data(data_paths):
    def __init__(self, file_name, no=1, method='', color='sep'):

        config = configparser.ConfigParser()
        config_file = 'config/' + file_name + '.ini'
        config.read(config_file)
        params = config['params']

        self.method         = method


        # read params
        self.path           = params.get('path')
        self.subdir         = params.get('subdir', None)
        self.color          = params.get('color', color)


        # sigma of gaussian_filter and threshold that is used in mask
        self.sigma          = params.getfloat('sigma', 5.)
        self.threshold      = params.getfloat('threshold', 1.3)
        # sigma that is used for halo substraction to improve mask
        self.halo_sig       = params.getfloat('halo_sig', None)

        # radius that is used in disk_filter and threshold to remove spots
        self.spots_radius   = params.getfloat('spots_radius', 10)
        self.thresh_spots   = params.getfloat('thresh_spots', 5.)


        self.lower_thresh   = params.getfloat('lower_thresh', None)

        # number of conncted areas that are interpreted as networks
        # '1' keeps only the largest area, None does not change the mask
        self.extract        = params.getint('extract', None)
        self.branch_thresh  = params.getint('branch_thresh', 50)


        self.analyse_flow   = params.getboolean('analyse_flow', False)
        self.symm_setup     = params.getboolean('symm_setup', False)


        # frame interval in seconds and scaling (pixel)
        self.frame_int      = params.getfloat('frame_int', 1.0)
        self.pixel_scaling  = params.getfloat('pixel_scaling', None)

        self.first          = params.getint('first', 1)
        self.last           = params.getint('last')

        # self.seed_positions = json.loads(params.get('seed_positions',
                                                    # '[[1000, 700]]'))

        image_prefix        = params.get('image_prefix')
        texas               = params.get('texas', '')
        green               = params.get('green', '')

        bf                  = params.get('bf', None)


        zeros               = params.getint('zeros', 3)
        super(data, self).__init__(image_prefix, no, zeros,
                                    texas, green,
                                    bf)
