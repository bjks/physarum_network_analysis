#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N animation
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -l h=!grotrian31
#$ -l h=!grotrian26
#$ -V
#$ -o /data.bpm/bksche/std/std_ani_out
#$ -e /data.bpm/bksche/std/std_ani_err
python3 /data.bpm/bksche/network_analysis/animation/animation_npz.py ${NAME} ${KEY} ${COL} ${NO}
