#! /bin/bash
# . ~/.bash_profile
# . ~/.bashrc
#$ -S /bin/bash
#$ -N phase
#$ -cwd
#$ -q grotrian.q
#$ -l h=!grotrian19
#$ -V
#$ -o /data.bpm/bksche/std/std_phase_out
#$ -e /data.bpm/bksche/std/std_phase_err
python3 /data.bpm/bksche/network_analysis/ratiometric/phase_hilbert.py ${NAME}
