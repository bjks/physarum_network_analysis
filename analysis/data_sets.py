import sys

class data_paths:
    def __init__(self, subdir=None):

        if subdir == None:
            dir = self.prefix
        else:
            dir = self.prefix + '/' + subdir

        if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
            if self.color != 'BF':
                self.file_raw1      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core1 + '_OG488.tif'
                self.file_raw2      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core2 + '_TexRe.tif'
            else:
                self.file_raw       = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core1  + '_Bright.tif'
            self.path_results   = '/data.bpm/bksche/results/' + self.prefix + '/'
            self.path_plots     = '/data.bpm/bksche/plots/' + self.prefix + '/'


        else:                               # local
            if self.color != 'BF':
                self.file_raw1      = '../image_data/'     + dir + '/' + self.core1 + '_OG488.tif'
                self.file_raw2      = '../image_data/'     + dir + '/' + self.core2 + '_TexRe.tif'
            else:
                self.file_raw       = '../image_data/' + dir + '/' + self.core1  + '_Bright.tif'
            self.path_results   = '../results/'    + self.prefix + '/'
            self.path_plots     = '../plots/'      + self.prefix + '/'


        self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
        self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

        self.file_plot      = self.path_plots + self.prefix + '_' + self.core1 + '_' + self.method + '_' + self.color
        self.file_plot_set  = self.path_plots + self.prefix + '_' + self.method + '_' + self.color


###################### 'my' data sets ####################
class data(data_paths):
    ### 1 (as eg in core1) corresponds to the actual dye (green)
    ### 2 corresponds to reference dye (red)
    def __init__(self, name, no=1, method='', color='gt'):
        if name == '2019_03_15':

            if color == 'gt' or color == 'green' or color == 'texas':
                self.core1           = 'TR_OG_t' + str(no).zfill(3)
                self.core2           = 'TR_OG_t' + str(no+1).zfill(3)
            elif color == 'tg':
                self.core1           = 'TR_OG_t' + str(no).zfill(3)
                self.core2           = 'TR_OG_t' + str(no).zfill(3)

            self.prefix         =   '2019-03-15'
            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 300
            self.thresh_spots   = 1.15
            self.spots_sig      = 20
            self.method         = method
            self.color          = color

            self.first, self.last = 1, 374
            self.seed_positions = [[500, 200], [250,1000]]

            super(data, self).__init__('TR_OG')


        if name == '2019_03_20_4':

            if color == 'gt' or color == 'green' or color == 'texas':
                self.core1           = 'TR_OG-04_t' + str(no).zfill(3)
                self.core2           = 'TR_OG-04_t' + str(no+1).zfill(3)
            elif color == 'tg':
                self.core1           = 'TR_OG-04_t' + str(no).zfill(3)
                self.core2           = 'TR_OG-04_t' + str(no).zfill(3)

            self.prefix         =   '2019-03-20'
            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 300
            self.thresh_spots   = 1.15
            self.spots_sig      = 20
            self.method         = method
            self.color          = color

            self.first, self.last = 1, 404
            self.seed_positions = [[500, 400]]

            super(data, self).__init__('TR_OG-04')

        if name == '2019_03_25':
            if color == 'gt' or color == 'green' or color == 'texas':
                self.core1       =   'TR_OG_t' + str(no).zfill(3)
                self.core2       =   'TR_OG_t' + str(no+1).zfill(3)
            elif color == 'tg':
                self.core1       =   'TR_OG_t' + str(no).zfill(3)
                self.core2       =   'TR_OG_t' + str(no).zfill(3)

            self.prefix         =   '2019-03-25'
            self.sigma          = 5
            self.threshold      = 1.3
            self.halo_sig       = None
            self.thresh_spots   = 1.15
            self.spots_sig      = 40
            self.method         = method
            self.color          = color

            self.first, self.last = 1, 310
            self.seed_positions = [[1500, 1000]]

            super(data, self).__init__('TR_OG')

        if name == '2019_03_29':
            if color == 'gt' or color == 'green' or color == 'texas':
                self.core1       =   'TR_OG_t' + str(no).zfill(4) ##### 4 !!
                self.core2       =   'TR_OG_t' + str(no+1).zfill(4)
            elif color == 'tg':
                self.core1       =   'TR_OG_t' + str(no).zfill(4) ##### 4 !!
                self.core2       =   'TR_OG_t' + str(no).zfill(4)

            self.prefix         =   '2019-03-29'
            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.15
            self.spots_sig      = 40.
            self.method         = method
            self.color          = color

            self.first, self.last = 1, 1141
            self.seed_positions = [[1500, 1500]]

            super(data, self).__init__('TR_OG')

        if name == '2019_04_11':
            if color == 'gt' or color == 'green' or color == 'texas':
                self.core1       =   'TR_OG-09_t' + str(no).zfill(3)
                self.core2       =   'TR_OG-09_t' + str(no+1).zfill(3)
            elif color == 'tg':
                self.core1       =   'TR_OG-09_t' + str(no).zfill(3)
                self.core2       =   'TR_OG-09_t' + str(no).zfill(3)
            elif color == 'BF':
                self.core1       =   'TR_OG-09_t' + str(no).zfill(3)

            self.prefix         =   '2019-04-11'
            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 30.
            self.method         = method
            self.color          = color

            self.first, self.last = 2, 849
            self.seed_positions = [[1000, 1000]]

            super(data, self).__init__('TR_OG-09')



####################################################################################
#################################  hard coded sets  ################################
####################################################################################
        if name == 'Mirna':
            self.prefix         = 'Mirna'
            if color == 'gt':
                self.core1           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no).zfill(4)+'c'
                self.core2           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no+1).zfill(4)+'c'
            elif color =='tg':
                self.core1           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no).zfill(4)+'c'
                self.core2           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no).zfill(4)+'c'

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = None
            self.thresh_spots   = 2.5
            self.spots_sig      = 10
            self.method         = method
            self.color          = color

            self.first, self.last = 1, 1005
            self.seed_positions = [[1000, 750], [750,1000], [1700, 1200], [600,750], [1500,900]]

            dir = self.prefix
            if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
                dir+= '/2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned'
                self.file_raw1      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core1    + '2.tif'
                self.file_raw2      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core2    + '1.tif'

                self.path_results   = '/data.bpm/bksche/results/' + self.prefix + '/'
                self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
                self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

                self.path_plots     = '/data.bpm/bksche/plots/' + self.prefix + '/'
                self.file_plot      = self.path_plots + self.prefix + '_' + self.core1 + '_' + self.method + '_' + self.color
                self.file_plot_set  = self.path_plots + self.prefix + '_' + self.method + '_' + self.color

            else:                               # local
                self.file_raw1      = 'image_data/'      + dir + '/' + self.core1   + '2.tif'
                self.file_raw2      = 'image_data/'      + dir + '/' + self.core2   + '1.tif'

                self.path_results   = '../results/'    + self.prefix + '/'
                self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
                self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

                self.path_plots     = '../plots/'      + self.prefix + '/'
                self.file_plot      = self.path_plots + self.prefix + '_' + self.core1 + '_' + self.method + '_' + self.color
                self.file_plot_set  = self.path_plots + self.prefix + '_' + self.method + '_' + self.color


        if name == 'SanDiego':
            self.prefix         = 'SanDiego'
            self.core           = 'img_'+str(no).zfill(8)+'0_shun'
            self.sigma          = 10
            self.threshold      = 1.1
            self.halo_sig       = None
            self.thresh_spots   = 1.5
            self.spots_sig      = 20
            self.method         = method
            self.first, self.last = 0, 240
            self.seed_positions = [[1000, 0]]

            dir = self.prefix

            if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
                self.file_raw1      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core + 'FITC_000' + '.tif'
                self.file_raw2      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core + 'TRITC_000' + '.tif'

                self.path_results   = '/data.bpm/bksche/results/' + self.prefix + '/'
                self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
                self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

                self.path_plots     = '/data.bpm/bksche/plots/' + self.prefix + '/'
                self.file_plot      = self.path_plots + self.prefix + '_' + self.core1 + '_' + self.method + '_' + self.color
                self.file_plot_set  = self.path_plots + self.prefix + '_' + self.method + '_' + self.color

            else:                               # local
                self.file_raw1      = 'image_data/'      + dir + '/' + self.core +  'FITC_000' + '.tif'
                self.file_raw2      = 'image_data/'      + dir + '/' + self.core +  'TRITC_000' + '.tif'

                self.path_results   = '../results/'    + self.prefix + '/'
                self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
                self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

                self.path_plots     = '../plots/'      + self.prefix + '/'
                self.file_plot      = self.path_plots + self.prefix + '_' + self.core1 + '_' + self.method + '_' + self.color
                self.file_plot_set  = self.path_plots + self.prefix + '_' + self.method + '_' + self.color
