#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N net_full
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_full_net_out
#$ -e /data.bpm/bksche/std/std_full_net_err
python3 /data.bpm/bksche/network_analysis/network/full_network.py ${NAME} ${COL} 20
