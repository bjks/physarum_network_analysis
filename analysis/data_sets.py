import sys
import os


#####################################################
################ Data Structure #####################
#####################################################

class data_paths:
    def __init__(self, green='_OG488', texas='_TexRe', bf='_Bright'):

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
        self.file_dat_set   = self.path_results + '/kymographs' + '_' + self.method + '_' + self.color

        self.file_plot      = self.path_plots + '_' + self.core1 + '_' + self.method + '_' + self.color
        self.file_plot_set  = self.path_plots + '_' + self.method + '_' + self.color



#####################################################
###################### data sets ####################
#####################################################


class data(data_paths):
    ### 1 (as eg in core1) corresponds to the actual dye (green)
    ### 2 corresponds to reference dye (red)
    def __init__(self, name, no=1, method='inter_mean', color='tg'):
        self.method         = method
        self.color          = color

        if name == 'Mirna':
            self.path         = 'Mirna'
            self.subdir       = None

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
            self.path         = 'SanDiego'
            self.subdir       = None

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
            self.path         = '2019-03-15'
            self.subdir       = 'TR_OG'

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

            super(data, self).__init__()


        if name == '2019_03_20_4':
            self.path         = '2019-03-20'
            self.subdir       = 'TR_OG-04'

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

            super(data, self).__init__()

        if name == '2019_03_25':
            self.path       = '2019-03-25'
            self.subdir     = 'TR_OG'

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

            super(data, self).__init__()

        if name == '2019_03_29':
            self.path       = '2019-03-29'
            self.subdir     = 'TR_OG'

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

            super(data, self).__init__()

        if name == '2019_04_11':
            self.path       = '2019-04-11'
            self.subdir     = 'TR_OG-09'

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

            super(data, self).__init__()

        if name == '2019_05_09_r':
            self.path       = '2019-05-09'
            self.subdir     = 'TR-CG1-4_AcquisitionBlock1_pt1'

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
            self.spots_sig      = 3.

            self.first, self.last = 1, 81
            self.seed_positions = [[750, 1500]]

            super(data, self).__init__(green='c2', texas='c1')


        if name == '2019_05_09_g':
            self.path       =   '2019-05-09'
            self.subdir     = 'TR-CG1-4_AcquisitionBlock2_pt2'

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
            self.spots_sig      = 3.

            self.first, self.last = 1, 2940
            self.seed_positions = [[750, 1500]]

            super(data, self).__init__(green='')

        if name == '2019_05_10_1_r':
            self.path          = '2019-05-10'
            self.subdir       = 'TR-CG1_AcquisitionBlock1_pt1'

            # TR-CG1_AcquisitionBlock1_pt1_t001c1.tif
            if color == 'gt':
                self.core1       =   'TR-CG1_AcquisitionBlock1_pt1_t' + str(no).zfill(3)
                self.core2       =   'TR-CG1_AcquisitionBlock1_pt1_t' + str(no+1).zfill(3)
            else:
                self.core1       =   'TR-CG1_AcquisitionBlock1_pt1_t' + str(no).zfill(3)
                self.core2       =   'TR-CG1_AcquisitionBlock1_pt1_t' + str(no).zfill(3)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 2.

            self.first, self.last = 1, 106
            self.seed_positions = [[800, 800]]

            super(data, self).__init__(green='c2', texas='c1')


        if name == '2019_05_10_2_r':
            self.path          = '2019-05-10'
            self.subdir       = 'TR-CG1_2_AcquisitionBlock1_pt1'

            # TR-CG1_AcquisitionBlock1_pt1_t001c1.tif
            if color == 'gt':
                self.core1       =   'TR-CG1_2_AcquisitionBlock1_pt1_t' + str(no).zfill(3)
                self.core2       =   'TR-CG1_2_AcquisitionBlock1_pt1_t' + str(no+1).zfill(3)
            else:
                self.core1       =   'TR-CG1_2_AcquisitionBlock1_pt1_t' + str(no).zfill(3)
                self.core2       =   'TR-CG1_2_AcquisitionBlock1_pt1_t' + str(no).zfill(3)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 2.

            self.first, self.last = 1, 106
            self.seed_positions = [[500, 250], [1000, 600]]

            super(data, self).__init__(green='c2', texas='c1')

        if name == '2019_05_10_2_g':
            self.path          = '2019-05-10'
            self.subdir       = 'TR-CG1_2_AcquisitionBlock2_pt2'

            # TR-CG1_AcquisitionBlock1_pt1_t001c1.tif
            if color == 'gt':
                self.core1       =   'TR-CG1_2_AcquisitionBlock2_pt2_t' + str(no).zfill(4)
                self.core2       =   'TR-CG1_2_AcquisitionBlock2_pt2_t' + str(no+1).zfill(4)
            else:
                self.core1       =   'TR-CG1_2_AcquisitionBlock2_pt2_t' + str(no).zfill(4)
                self.core2       =   'TR-CG1_2_AcquisitionBlock2_pt2_t' + str(no).zfill(4)

            self.sigma          = 5
            self.threshold      = 1.2
            self.halo_sig       = 200.
            self.thresh_spots   = 1.3
            self.spots_sig      = 2.

            self.first, self.last = 1, 1136
            self.seed_positions = [[500, 250], [1000, 600]]

            super(data, self).__init__(green='')
