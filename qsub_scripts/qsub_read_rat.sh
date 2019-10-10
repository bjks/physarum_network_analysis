#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N read_rat
#$ -cwd
#$ -q grotrian.q
#$ -q mvapich2-grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_read_out
#$ -e /data.bpm/bksche/std/std_read_err
python3 /data.bpm/bksche/network_analysis/read_ratiometric.py ${NAME} ${START} ${STOP} ${STEP}
