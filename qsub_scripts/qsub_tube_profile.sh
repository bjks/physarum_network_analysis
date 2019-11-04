#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N tube_profile
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_tube_out
#$ -e /data.bpm/bksche/std/std_tube_err
python3 /data.bpm/bksche/network_analysis/tube_profile.py ${NAME} ${NO}
