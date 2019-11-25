#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N rat_full
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_full_out
#$ -e /data.bpm/bksche/std/std_full_err
python3 /data.bpm/bksche/network_analysis/ratiometric/full_ratiometric.py ${NAME} 20
