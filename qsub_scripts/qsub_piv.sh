#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N piv
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_piv_out
#$ -e /data.bpm/bksche/std/std_piv_err
python3 /data.bpm/bksche/network_analysis/piv/piv.py ${NAME} ${COL} ${START} ${END}
