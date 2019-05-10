import os
import sys
from analysis.data_sets import *

to_qsub = ' ' + '/data.bpm/bksche/network_analysis/qsub_scripts/' # space at beginning


#####
def command_ratiometric(name, a, b):
    return 'qsub -v NAME=' + name + ',START=' + str(a) + ',END=' + str(b) + to_qsub + 'qsub_ratiometric.sh'

def command_skeleton(name, no):
    return 'qsub -v NAME=' + name + ',NO=' + str(no) + to_qsub + 'qsub_skeleton.sh'

#####
def command_network(name, a, b, col):
    return 'qsub -v NAME=' + name + ',COL=' + str(col) + ',START=' + str(a) + ',END=' + str(b) + to_qsub + 'qsub_network.sh'

def command_skeleton_network(name, no, col):
    return 'qsub -v NAME=' + name + ',COL=' + str(col) + ',NO=' + str(no) + to_qsub + 'qsub_skeleton_network.sh'


#####
def command_correlation(name, no):
    return 'qsub -v NAME=' + name + ',NO=' + str(no) + to_qsub + 'qsub_correlation.sh'

def command_tube_profile(name, no):
    return 'qsub -v NAME=' + name + ',NO=' + str(no) + to_qsub + 'qsub_tube_profile.sh'

def command_animation(name, no):
    return 'qsub -v NAME=' + name + ',NO=' + str(no) + to_qsub +  'qsub_animation.sh'

def run_ratiometric():
    name        = os.sys.argv[2]
    no_slices   = int(os.sys.argv[3])
    set         = data(name)
    first, last = data(name).first, data(name).last

    ######################
    ######################

    last+=1
    total = last-first

    set_size = int( total/no_slices )
    last_sub = set_size * no_slices + first

    ### distr sets to slices
    for i in range(0, no_slices):
        start   = set_size * i +first
        end     = set_size * (i+1) + first
        print(command_ratiometric(name, start, end))
        if sys.platform.startswith('linux'):
            os.system(command_ratiometric(name, start, end))

    print('remaining: ')
    ### if some remain...
    while last_sub < last:
        print(command_ratiometric(name, last_sub, last_sub+1))
        if sys.platform.startswith('linux'):
            os.system(command_ratiometric(name, last_sub, last_sub+1))
        last_sub +=1



def run_skeleton():
    name                = os.sys.argv[2]
    set                 = data(name)
    no_seed_positions   = len(data(name).seed_positions)

    for i in range(no_seed_positions):
        print('Seed index: ', i)
        if sys.platform.startswith('linux'):
            os.system(command_skeleton(name, i))

def run_skeleton_network():
    name                = os.sys.argv[2]
    col                 = os.sys.argv[3]

    set                 = data(name)
    no_seed_positions   = len(data(name).seed_positions)

    for i in range(no_seed_positions):
        print('Network: seed index: ', i)
        if sys.platform.startswith('linux'):
            os.system(command_skeleton_network(name, i, col))


def run_correlation():
    name                = os.sys.argv[2]
    set                 = data(name)

    print(command_correlation(name, int(os.sys.argv[3])))
    if sys.platform.startswith('linux'):
        os.system(command_correlation(name, int(os.sys.argv[3])))


def run_tube_profile():
    name                = os.sys.argv[2]
    set                 = data(name)

    print(command_tube_profile(name, int(os.sys.argv[3])))
    if sys.platform.startswith('linux'):
        os.system(command_tube_profile(name, int(os.sys.argv[3])))

def run_animation():
    name                = os.sys.argv[2]
    set                 = data(name)

    print(command_animation(name, int(os.sys.argv[3])))
    if sys.platform.startswith('linux'):
        os.system(command_animation(name, int(os.sys.argv[3])))



def run_network():
    name        = os.sys.argv[2]
    no_slices   = int(os.sys.argv[3])
    col         = os.sys.argv[4]
    set         = data(name)
    first, last = data(name).first, data(name).last

    ######################
    ######################

    last+=1
    total = last-first

    set_size = int( total/no_slices )
    last_sub = set_size * no_slices + first

    ### distr sets to slices
    for i in range(0, no_slices):
        start   = set_size * i +first
        end     = set_size * (i+1) + first
        print(command_network(name, start, end, col))
        if sys.platform.startswith('linux'):
            os.system(command_network(name, start, end, col))

    print('remaining: ')
    ### if some remain...
    while last_sub < last:
        print(command_network(name, last_sub, last_sub+1, col))
        if sys.platform.startswith('linux'):
            os.system(command_network(name, last_sub, last_sub+1, col))
        last_sub +=1

# if 'ratiometric'  : ...
# if 'skeleton'     : ...
# if 'correlation'  : ...


script      = os.sys.argv[1]

if script.startswith('rat'):
    run_ratiometric()

elif script.startswith('ske'):
    run_skeleton()

elif script.startswith('cor'):
    run_correlation()

elif script.startswith('net'):
    run_network()

elif script.startswith('tube'):
    run_tube_profile()

elif script.startswith('ani'):
    run_animation()

elif script.startswith('s_n'):
    run_skeleton_network()
