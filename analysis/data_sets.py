import sys
import os


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

        if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
            dir_raw             = '/data.bpm/zeiss/BjoernZ/' + dir + '/'
            dir_results         = '/data.bpm/bksche/results/'
            dir_plots           = '/data.bpm/bksche/plots/'

        else:                               # local
            dir_raw             = '../image_data/'     + dir + '/'
            dir_results         = '../results/'
            dir_plots           = '../plots/'


        self.file_raw       = dir_raw + self.core1 + bf + '.tif'
        self.file_raw1      = dir_raw + self.core1 + green + '.tif'
        self.file_raw2      = dir_raw + self.core2 + texas + '.tif'


        for d in [dir_results, dir_plots]:
            if not os.path.exists(d + self.path):
                os.mkdir(d + self.path)
            if self.subdir != None:
                if not os.path.exists(d + self.path + '/' +self.subdir):
                    os.mkdir(d + self.path + '/' +self.subdir)


        self.path_results   = '../results/'    + dir + '/'
        self.path_plots     = '../plots/'      + dir + '/'

        self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
        self.file_dat_set   = self.path_results + 'kymographs' + '_' + self.method + '_' + self.color

        self.file_plot      = self.path_plots + self.core1 + '_' + self.method + '_' + self.color
        self.file_plot_set  = self.path_plots + self.method + '_' + self.color



#####################################################
###################### data sets ####################
#####################################################
import configparser
import json

class data(data_paths):
    def __init__(self, file_name, no=1, method='inter_mean', color='tg'):

        config = configparser.ConfigParser()
        config_file = 'config/' + file_name + '.ini'
        config.read(config_file)
        params = config['params']

        self.method         = method
        self.color          = color

        self.path           = params.get('path')
        self.subdir         = params.get('subdir', None)
        self.sigma          = float(params.get('sigma', 5.))
        self.threshold      = float(params.get('thresh_spots', 1.3))
        self.halo_sig       = float(params.get('halo_sig', 0))
        self.thresh_spots   = float(params.get('thresh_spots', 5.))
        self.spots_sig      = float(params.get('spots_sig', 10))

        self.first          = int(params.get('first', 1))
        self.last           = int(params.get('last'))

        self.seed_positions = json.loads(params.get('seed_positions'))

        image_prefix        = params.get('image_prefix')
        texas               = params.get('texas', '')
        green               = params.get('green', '')
        bf                  = params.get('bf', '')
        zeros               = int(params.get('zeros', 3))
        super(data, self).__init__(image_prefix, no, zeros, texas, green, bf)
