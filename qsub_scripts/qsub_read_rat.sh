#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N read_rat
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -l h=!grotrian31
#$ -l h=!grotrian26
#$ -V
#$ -o /data.bpm/bksche/std/std_read_out
#$ -e /data.bpm/bksche/std/std_read_err
python3 /data.bpm/bksche/network_analysis/ratiometric/read_ratiometric.py ${NAME} ${START} ${STOP} ${STEP}
