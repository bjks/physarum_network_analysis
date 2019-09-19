#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N skeleton
#$ -cwd
#$ -q mvapich2-grotrian.q
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_skel_out
#$ -e /data.bpm/bksche/std/std_skel_err
python3 /data.bpm/bksche/network_analysis/skeleton.py ${NAME}
