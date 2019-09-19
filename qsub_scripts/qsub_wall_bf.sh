#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N wall_bf
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -q mvapich2-grotrian.q
#$ -V
#$ -o /data.bpm/bksche/std/std_wall_out
#$ -e /data.bpm/bksche/std/std_wall_err
python3 /data.bpm/bksche/network_analysis/wall_bf.py 2019-08-13_phyta_over2 1
