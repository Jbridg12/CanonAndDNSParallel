#PBS -S /bin/bash
#PBS -A ACF-UTK0175
#PBS -o $HOME/$PBS_JOBID-out.txt
#PBS -l partition=general
#PBS -l nodes=24:ppn=22
#PBS -l walltime=1:00:00

mpicc -lm $HOME/canon.c
mpirun -n 1 $HOME/a.out 1024 1
mpirun -n 16 $HOME/a.out 1024 16
mpirun -n 64 $HOME/a.out 1024 64
mpirun -n 256 $HOME/a.out 1024 256