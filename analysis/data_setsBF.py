import sys

class data_pathsBF:
    def __init__(self, subdir=None):

        if subdir == None:
            dir = self.prefix
        else:
            dir = self.prefix + '/' + subdir

        if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
            self.file_raw      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core + '_Bright.tif'

        else:                               # local
            self.file_raw      = 'image_data/'     + dir + '/' + self.core + '_Bright.tif'

        self.path_results   = 'results/'    + self.prefix + '/'
        self.file_dat       = 'results/'    + self.prefix + '/' + self.core + '_' + self.method
        self.file_dat_set   = 'results/'    + self.prefix + '/kymographs'   + '_' + self.method

        self.path_plots     = 'plots/'      + self.prefix + '/'
        self.file_plot      = 'plots/'      + self.prefix + '/' + self.prefix + '_' + self.core + '_' + self.method
        self.file_plot_set  = 'plots/'      + self.prefix + '/' + self.prefix + '_' + self.method


###################### 'my' data sets ####################
class dataBF(data_pathsBF):
    def __init__(self, name, no=1, method='', order='gt'):
        if name == '2019_04_11':
            # TR_OG-09_t001_Bright.tif

            self.core           = 'TR_OG-09_t' + str(no).zfill(3)

            self.prefix         =   '2019-04-11'
            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 300
            self.thresh_spots   = 1.15
            self.spots_sig      = 20
            self.method         = method

            self.first, self.last = 1, 858
            self.seed_positions = [[1000, 1000]]

            super(dataBF, self).__init__('TR_OG-09')
