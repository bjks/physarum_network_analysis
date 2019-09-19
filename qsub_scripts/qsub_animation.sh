#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N animation
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -q mvapich2-grotrian.q
#$ -V
#$ -o /data.bpm/bksche/std/std_ani_out
#$ -e /data.bpm/bksche/std/std_ani_err
python3 /data.bpm/bksche/network_analysis/animation_conc.py ${NAME} ${KEY} ${COL} ${NO}
