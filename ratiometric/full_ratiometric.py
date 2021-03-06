import os
import time
import subprocess
import datetime

import sys
sys.path.append("..")

from analysis.data_sets import *

############################################
############################################
NAME        = os.sys.argv[1]
NO_SLICES   = int(os.sys.argv[2])
COLOR       = 'sep'
MAX_FRAME   = 200
SKIP        = False
WAIT_FORCZI = False      # wait for running czi converters
############################################
############################################

to_qsub = ' ' + '/data.bpm/bksche/network_analysis/qsub_scripts/' # space at beginning

############################################

def sub_command(com):
    print(com)
    if sys.platform.startswith('linux'):
        os.system(com)
    else:
        print("Dryrun!")


def command_ratiometric(a, b):
    return 'qsub -v NAME=' + NAME + ',START=' + str(a) + ',END=' + str(b) \
            + to_qsub + 'qsub_ratiometric.sh'

def run_ratiometric():
    name        = NAME
    no_slices   = NO_SLICES
    first, last = data(NAME).first, data(NAME).last

    last += 1
    total = last-first

    set_size = int( total/no_slices )
    last_sub = set_size * no_slices + first

    ### distr sets to slices
    for i in range(0, no_slices):
        start   = set_size * i +first
        end     = set_size * (i+1) + first
        com = command_ratiometric(start, end)
        sub_command(com)

    print('remaining: ')
    ### if some remain...
    while last_sub < last:
        com = command_ratiometric(last_sub, last_sub+1)
        sub_command(com)
        last_sub += 1

############################################

def run_skeleton():
    com  =  'qsub -v NAME=' + NAME + ',COL=' + COLOR + to_qsub + 'qsub_skeleton.sh'
    sub_command(com)

def run_snr():
    com  =  'qsub -v NAME=' + NAME + to_qsub + 'qsub_snr.sh'
    sub_command(com)

def run_animation(key):
    first, last = data(NAME).first, data(NAME).last
    step = int((last - first)/MAX_FRAME)
    com = 'qsub -v NAME=' + NAME + \
                  ',KEY=' + key + \
                  ',COL=' + COLOR + \
                  ',NO='  + str(step) + to_qsub +  'qsub_animation.sh'
    sub_command(com)


def run_tube_profile():
    com  = 'qsub -v NAME=' + NAME + \
                    ',NO=' + str(2) + to_qsub + 'qsub_tube_profile.sh'
    sub_command(com)


def run_read():
    first, last = data(NAME).first, data(NAME).last
    com = 'qsub -v NAME=' + NAME + \
                  ',START=' + str(first) + \
                  ',STOP=' + str(last) + \
                  ',STEP='  + str(50) + to_qsub +  'qsub_read_rat.sh'
    sub_command(com)

############################################
def run_phase():
    com  =  'qsub -v NAME=' + NAME + to_qsub + 'qsub_phase.sh'
    sub_command(com)



###################################################
def check_and_wait(sec, script):
    while True and sys.platform.startswith('linux'):
        # print("Checking qstat...")
        qstat = subprocess.check_output("qstat", universal_newlines=True)
        # print(len(qstat))
        if script not in qstat:
            break
        time.sleep(sec)

####################################################################
####################################################################
print("Start analysis of ", NAME, " at ", datetime.datetime.now())

if  not SKIP:
    if WAIT_FORCZI:
        check_and_wait(100, 'czi_')
    run_ratiometric()
    check_and_wait(100, 'ratiomet')


run_skeleton()
run_snr()
run_tube_profile()
run_read()
run_animation('skeleton')
run_animation('raw')


check_and_wait(100, 'skelet')
run_phase()

print("Phase analysis of ", NAME, " submitted at ", datetime.datetime.now())
