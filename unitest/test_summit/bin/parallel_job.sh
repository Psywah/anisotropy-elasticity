#!/bin/bash

#SBATCH --job-name carotidIso
#SBATCH --qos normal
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 24
#SBATCH --time 2:00:00
#SBATCH --output carotid.out
#SBATCH --partition shas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=psywah@gmail.com

## SBATCH --ntasks 48

# the slurm module provides the srun command

cd /scratch/summit/shgo7817/anisotropy-elasticity/unitest/test_summit/bin
source my-config

#mpiexec -n 4 ./CarotidIso.ex
mpirun -n 48 ./CarotidIso.ex
#echo Running on $(hostname --fqdn):  'Hello, world!'
