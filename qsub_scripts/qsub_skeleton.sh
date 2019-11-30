#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N skeleton
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -l h=!grotrian31
#$ -l h=!grotrian26
#$ -V
#$ -o /data.bpm/bksche/std/std_skel_out
#$ -e /data.bpm/bksche/std/std_skel_err
python3 /data.bpm/bksche/network_analysis/skeleton/skeleton.py ${NAME} ${COL}
