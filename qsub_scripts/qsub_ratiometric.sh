#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N ratiometric
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -l h=!grotrian31
#$ -l h=!grotrian26
#$ -V
#$ -o /data.bpm/bksche/std/std_rat_out
#$ -e /data.bpm/bksche/std/std_rat_err
python3 /data.bpm/bksche/network_analysis/ratiometric/ratiometric.py ${NAME} ${START} ${END}
