#PBS -S /bin/bash
#PBS -A ACF-UTK0175
#PBS -o $HOME/$PBS_JOBID-out.txt
#PBS -l partition=general
#PBS -l nodes=24:ppn=22
#PBS -l walltime=1:00:00

mpicc -lm $HOME/dns.c
mpirun -n 1 $HOME/a.out 1024 1
mpirun -n 8 $HOME/a.out 1024 8
mpirun -n 64 $HOME/a.out 1024 64
mpirun -n 512 $HOME/a.out 1024 512