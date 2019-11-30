#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N phase_net
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -l h=!grotrian31
#$ -l h=!grotrian26
#$ -V
#$ -o /data.bpm/bksche/std/std_phase_net_out
#$ -e /data.bpm/bksche/std/std_phase_net_err
python3 /data.bpm/bksche/network_analysis/network/phase_hilbert_net.py ${NAME} ${COL}
