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
    """
    Takes string parameters and composes filenames of raw images and output
    files
    Note that this class does not know about config files
    """
    def __init__(self, image_prefix, no, zeros, texas, green, bf):

        if self.subdir == None:
            dir = self.path
        else:
            dir = self.path + '/' + self.subdir

        if sys.platform.startswith('linux'): # i.e. on the (ubuntu) cluster
            dir_raw             = '/data.bpm/zeiss/BjoernZ/' + dir + '/'
            dir_results         = '/data.bpm/bksche/results/'
            dir_plots           = '/data.bpm/bksche/plots/'
            dir_logs            = '/data.bpm/bksche/logs/'

        else:                               # local
            dir_raw             = '/Users/bjoern/image_analysis/image_data/' + dir + '/'
            dir_results         = '/Users/bjoern/image_analysis/results/'
            dir_plots           = '/Users/bjoern/image_analysis/plots/'
            dir_logs            = '/Users/bjoern/image_analysis/logs/'


        for d in [dir_results, dir_plots]:
            mk_mising_dir(d + self.path)
            if self.subdir != '':
                mk_mising_dir(d + self.path + '/' +self.subdir)

        core = image_prefix + str(no).zfill(zeros)

        if bf != None:
            bf = [x.strip() for x in bf.split()]
            if np.size(bf)>1:
                self.file_raw = [dir_raw + core + b + '.tif' for b in bf]
            elif np.size(bf)==1:
                self.file_raw = dir_raw + core + bf[0] + '.tif'
            else:
                self.file_raw = dir_raw + core + '.tif'

        if green != None:
            self.file_raw1      = dir_raw + core + green + '.tif'
        if texas != None:
            self.file_raw2      = dir_raw + core + texas + '.tif'


        self.path_results   = dir_results + dir + '/'
        self.path_plots     = dir_plots + dir + '/'

        full_file = core + '_' + self.method + '_' + self.color
        set_file  = 'set_' + self.method + '_' + self.color

        if self.lower_thresh != None:
            full_file += '_back_corr'
            set_file  += '_back_corr'

        self.file_dat       = self.path_results + full_file
        self.file_dat_flow  = self.path_results + full_file + '_flow'

        self.file_dat_set   = self.path_results + set_file

        self.file_plot              = mk_mising_dir(self.path_plots + 'samples/') + full_file
        self.file_plot_tube_profile = mk_mising_dir(self.path_plots + 'tube_profile/') + set_file
        self.file_plot_set          = self.path_plots + set_file

        self.file_log = dir_logs + self.keyword + '.log'


#####################################################
###################### data sets ####################
#####################################################
import configparser
import json

class data(data_paths):
    """
    Class interpreting config files given as keyword and read as
    '../config/' + keyword + '.ini' (relative path!)
    calls super class constructor data_paths which handles filenames and
    data structure
    """
    def __init__(self, keyword, no=1, method='', color='sep'):

        config = configparser.ConfigParser()
        config_file = '../config/' + keyword + '.ini'
        config.read(config_file)
        params = config['params']

        self.keyword        = keyword
        self.method         = method


        # directories. subdir is optional
        self.path           = params.get('path')
        self.subdir         = params.get('subdir', '')
        self.color          = params.get('color', color)


        # sigma of gaussian_filter and threshold that is used in mask
        self.sigma          = params.getfloat('sigma', 5.)
        self.threshold      = params.getfloat('threshold', 1.3)
        # sigma that is used for halo substraction to improve mask, optional
        self.halo_sig       = params.getfloat('halo_sig', None)

        # radius that is used in disk_filter and threshold to remove spots
        # if one of these params is not given, remove_spots is scipped
        self.spots_radius   = params.getfloat('spots_radius', None)
        self.thresh_spots   = params.getfloat('thresh_spots', None)

        # lower threshold to get only background, used for background_correction
        self.lower_thresh   = params.getfloat('lower_thresh', None)

        # number of conncted areas that are interpreted as networks
        # '1' keeps only the largest area, None does not change the mask
        self.extract        = params.getint('extract', 1)
        self.branch_thresh  = params.getint('branch_thresh', 10)

        # if analyse flow, flow_analysis is called, bf frames as channels is
        # assumed
        self.analyse_flow   = params.getboolean('analyse_flow', False)

        # if green and texas channel are images in a symmetric way,
        # interpolation is used to shift signals accordingly to account for
        # time shift between imaging
        # assumed order per frame: 1. texas 2. green
        self.symm_setup     = params.getboolean('symm_setup', False)

        # frame interval in seconds and scaling (um(!) per pixel)
        self.frame_int      = params.getfloat('frame_int', 1.0)
        self.pixel_scaling  = params.getfloat('pixel_scaling', None)

        self.first          = params.getint('first', 1)
        self.last           = params.getint('last')

        # times (in frame number) and positions (in pixels) can be defined
        # to focus only on these regions in kymograph_analysis
        self.times = tuple(int(t) if t != 'None' else None
                        for t in params.get('times', 'None None').split())

        self.positions = tuple(int(t) if t != 'None' else None
                        for t in params.get('positions', 'None None').split())


        # channel naming, image prefix (stuff before frame number)
        image_prefix        = params.get('image_prefix', self.subdir + '_t')
        texas               = params.get('texas', None)
        green               = params.get('green', None)

        bf                  = params.get('bf', None)

        # digits in image names, note that this is not given by total number of
        # frames as analysis might only handle a subset of available images
        zeros               = params.getint('zeros', 3)

        super(data, self).__init__(image_prefix, no, zeros,
                                    texas, green,
                                    bf)
