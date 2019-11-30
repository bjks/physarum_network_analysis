#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N full_flow
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_full_piv_out
#$ -e /data.bpm/bksche/std/std_full_piv_err
python3 /data.bpm/bksche/network_analysis/piv/full_piv.py ${NAME} bf 20
