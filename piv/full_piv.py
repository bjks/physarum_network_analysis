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
COLOR       = os.sys.argv[2]
NO_SLICES   = int(os.sys.argv[3])
MAX_FRAME   = 200
SKIP        = False
IGNORE_CONV = True      # ignore running czi converters
############################################
############################################

to_qsub = ' ' + '/data.bpm/bksche/network_analysis/qsub_scripts/' # space at beginning

############################################

def command_piv(a, b):
    return 'qsub -v NAME=' + NAME + ',COL' + COLOR + ',START=' + str(a) + \
            ',END=' + str(b) + to_qsub + 'qsub_piv.sh'

def sub_command(com):
    print(com)
    if sys.platform.startswith('linux'):
        os.system(com)
    else:
        print("Dryrun!")

def run_piv():
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
        com = command_piv(start, end)
        sub_command(com)

    print('remaining: ')
    ### if some remain...
    while last_sub < last:
        com = command_piv(last_sub, last_sub+1)
        sub_command(com)
        last_sub += 1


def run_animation(key):
    first, last = data(NAME).first, data(NAME).last
    step = int((last - first)/MAX_FRAME)
    com = 'qsub -v NAME=' + NAME + \
                  ',KEY=' + key + \
                  ',COL=' + COLOR + \
                  ',NO='  + str(step) + to_qsub +  'qsub_animation.sh'
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

if SKIP:
    pass
else:
    if not IGNORE_CONV:
        check_and_wait(100, 'czi_')
    run_piv()
    check_and_wait(100, 'piv')

run_animation('flow')
run_animation('mask')
