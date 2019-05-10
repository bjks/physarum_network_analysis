import sys


#####################################################
################ Data Structure #####################
#####################################################

class data_paths:
    def __init__(self, subdir=None, green='_OG488', texas='_TexRe', bf='_Bright'):

        if subdir == None:
            dir = self.prefix
        else:
            dir = self.prefix + '/' + subdir

        if sys.platform.startswith('linux'): # i.e. you are on the (ubuntu) cluster
            self.file_raw1      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core1 + green + '.tif'
            self.file_raw2      = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core2 + texas + '.tif'
            self.file_raw       = '/data.bpm/zeiss/BjoernZ/' + dir + '/' + self.core1 + bf + '.tif'

            self.path_results   = '/data.bpm/bksche/results/'   + self.prefix + '/'
            self.path_plots     = '/data.bpm/bksche/plots/'     + self.prefix + '/'


        else:                               # local
            self.file_raw1      = '../image_data/'     + dir + '/' + self.core1 + green + '.tif'
            self.file_raw2      = '../image_data/'     + dir + '/' + self.core2 + texas + '.tif'
            self.file_raw       = '../image_data/'     + dir + '/' + self.core1 + bf + '.tif'

            self.path_results   = '../results/'    + self.prefix + '/'
            self.path_plots     = '../plots/'      + self.prefix + '/'


        self.file_dat       = self.path_results + self.core1    + '_' + self.method + '_' + self.color
        self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

        self.file_plot      = self.path_plots + self.prefix + '_' + self.core1 + '_' + self.method + '_' + self.color
        self.file_plot_set  = self.path_plots + self.prefix + '_' + self.method + '_' + self.color



#####################################################
###################### data sets ####################
#####################################################


class data(data_paths):
    ### 1 (as eg in core1) corresponds to the actual dye (green)
    ### 2 corresponds to reference dye (red)
    def __init__(self, name, no=1, method='', color='gt'):
        self.method         = method
        self.color          = color

        if name == 'Mirna':
            self.prefix         = 'Mirna'

            if color == 'gt':
                self.core1           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no).zfill(4)
                self.core2           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no+1).zfill(4)
            else:
                self.core1           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no).zfill(4)
                self.core2           = '2016-09-22_5uL_3.3mM_TR_10uL_0.1mMORG_pruned_t'+str(no).zfill(4)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = None
            self.thresh_spots   = 2.5
            self.spots_sig      = 10

            self.first, self.last = 1, 1005
            self.seed_positions = [[1000, 750], [750,1000], [1700, 1200], [600,750], [1500,900]]
            super(data, self).__init__(green='c2', texas='c1')


        if name == 'SanDiego':
            self.prefix         = 'SanDiego'

            if color == 'gt':
                self.core1           = 'img_'+str(no).zfill(8)+'0_shun'
                self.core2           = 'img_'+str(no).zfill(8)+'0_shun'
            else:
                self.core1           = 'img_'+str(no+1).zfill(8)+'0_shun'
                self.core2           = 'img_'+str(no).zfill(8)+'0_shun'

            self.sigma          = 10
            self.threshold      = 1.1
            self.halo_sig       = None
            self.thresh_spots   = 1.5
            self.spots_sig      = 20

            self.first, self.last = 0, 240
            self.seed_positions = [[1000, 0]]
            super(data, self).__init__(green='FITC_000', texas='TRITC_000')


        if name == '2019_03_15':
            self.prefix         =   '2019-03-15'

            if color == 'gt':
                self.core1           = 'TR_OG_t' + str(no).zfill(3)
                self.core2           = 'TR_OG_t' + str(no+1).zfill(3)
            else:
                self.core1           = 'TR_OG_t' + str(no).zfill(3)
                self.core2           = 'TR_OG_t' + str(no).zfill(3)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 300
            self.thresh_spots   = 1.15
            self.spots_sig      = 20

            self.first, self.last = 1, 374
            self.seed_positions = [[500, 200], [250,1000]]

            super(data, self).__init__('TR_OG')


        if name == '2019_03_20_4':
            self.prefix         =   '2019-03-20'

            if color == 'gt':
                self.core1           = 'TR_OG-04_t' + str(no).zfill(3)
                self.core2           = 'TR_OG-04_t' + str(no+1).zfill(3)
            else:
                self.core1           = 'TR_OG-04_t' + str(no).zfill(3)
                self.core2           = 'TR_OG-04_t' + str(no).zfill(3)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 300
            self.thresh_spots   = 1.15
            self.spots_sig      = 20

            self.first, self.last = 1, 404
            self.seed_positions = [[500, 400]]

            super(data, self).__init__('TR_OG-04')

        if name == '2019_03_25':
            self.prefix         =   '2019-03-25'

            if color == 'gt':
                self.core1       =   'TR_OG_t' + str(no).zfill(3)
                self.core2       =   'TR_OG_t' + str(no+1).zfill(3)
            else:
                self.core1       =   'TR_OG_t' + str(no).zfill(3)
                self.core2       =   'TR_OG_t' + str(no).zfill(3)

            self.sigma          = 5
            self.threshold      = 1.3
            self.halo_sig       = None
            self.thresh_spots   = 1.15
            self.spots_sig      = 40

            self.first, self.last = 1, 310
            self.seed_positions = [[1500, 1000]]

            super(data, self).__init__('TR_OG')

        if name == '2019_03_29':
            self.prefix         =   '2019-03-29'

            if color == 'gt':
                self.core1       =   'TR_OG_t' + str(no).zfill(4) ##### 4 !!
                self.core2       =   'TR_OG_t' + str(no+1).zfill(4)
            else:
                self.core1       =   'TR_OG_t' + str(no).zfill(4) ##### 4 !!
                self.core2       =   'TR_OG_t' + str(no).zfill(4)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.15
            self.spots_sig      = 40.

            self.first, self.last = 1, 1141
            self.seed_positions = [[1500, 1500]]

            super(data, self).__init__('TR_OG')

        if name == '2019_04_11':
            self.prefix         =   '2019-04-11'

            if color == 'gt':
                self.core1       =   'TR_OG-09_t' + str(no).zfill(3)
                self.core2       =   'TR_OG-09_t' + str(no+1).zfill(3)
            else:
                self.core1       =   'TR_OG-09_t' + str(no).zfill(3)
                self.core2       =   'TR_OG-09_t' + str(no).zfill(3)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 30.

            self.first, self.last = 2, 849
            self.seed_positions = [[1000, 1000]]

            super(data, self).__init__('TR_OG-09')

        if name == '2019_05_09_1':
            self.prefix         =   '2019-05-09'

            # TR-CG1-4_AcquisitionBlock1_pt1_t82c2.tif
            if color == 'gt':
                self.core1       =   'TR-CG1-4_AcquisitionBlock1_pt1_t' + str(no).zfill(2)
                self.core2       =   'TR-CG1-4_AcquisitionBlock1_pt1_t' + str(no+1).zfill(2)
            else:
                self.core1       =   'TR-CG1-4_AcquisitionBlock1_pt1_t' + str(no).zfill(2)
                self.core2       =   'TR-CG1-4_AcquisitionBlock1_pt1_t' + str(no).zfill(2)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 30.

            self.first, self.last = 1, 81
            self.seed_positions = [[750, 1500]]

            super(data, self).__init__(subdir='TR-CG1-4_AcquisitionBlock1_pt1', green='c2', texas='c1')


        if name == '2019_05_09_2':
            self.prefix         =   '2019-05-09'

            # TR-CG1-4_AcquisitionBlock1_pt1_t82c2.tif
            if color == 'gt':
                self.core1       =   'TR-CG1-4_AcquisitionBlock2_pt2_t' + str(no).zfill(4)
                self.core2       =   'TR-CG1-4_AcquisitionBlock2_pt2_t' + str(no+1).zfill(4)
            else:
                self.core1       =   'TR-CG1-4_AcquisitionBlock2_pt2_t' + str(no).zfill(4)
                self.core2       =   'TR-CG1-4_AcquisitionBlock2_pt2_t' + str(no).zfill(4)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 30.

            self.first, self.last = 1, 2940
            self.seed_positions = [[750, 1500]]

            super(data, self).__init__(subdir='TR-CG1-4_AcquisitionBlock2_pt2', green='')
