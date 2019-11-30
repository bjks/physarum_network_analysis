#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N network
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -l h=!grotrian31
#$ -l h=!grotrian26
#$ -V
#$ -o /data.bpm/bksche/std/std_net_out
#$ -e /data.bpm/bksche/std/std_net_err
python3 /data.bpm/bksche/network_analysis/network/network.py ${NAME} ${COL} ${START} ${END}
