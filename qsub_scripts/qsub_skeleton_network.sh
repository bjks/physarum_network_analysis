#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N skeleton_net
#$ -cwd
#$ -q mvapich2-grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_skel_n_out
#$ -e /data.bpm/bksche/std/std_skel_n_err
python3 /data.bpm/bksche/network_analysis/skeleton_network.py ${NAME} ${COL} ${NO}
