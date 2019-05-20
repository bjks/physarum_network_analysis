#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N corr
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std_out
#$ -e /data.bpm/bksche/std_err
python3 /data.bpm/bksche/network_analysis/correlation.py ${NAME} ${COL}
