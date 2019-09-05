#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N rat_full
#$ -cwd
#$ -q mvapich2-grotrian.q
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_full_out
#$ -e /data.bpm/bksche/std/std_full_err
python3 /data.bpm/bksche/network_analysis/full_analysis.py 2019-08-08 50
