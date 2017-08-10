#!/bin/bash

#SBATCH --job-name fineMeshanalyticJac
#SBATCH --qos normal
#SBATCH --nodes 10
#SBATCH --ntasks-per-node 24
#SBATCH --time 2:00:00
#SBATCH --output analyticJac-finemesh.out
#SBATCH --partition shas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=psywah@gmail.com

## SBATCH --ntasks 240

# the slurm module provides the srun command

cd /scratch/summit/shgo7817/local-test/analyticJacRes/bin
source my-config

#mpiexec -n 4 ./CarotidIso.ex
mpirun -n 240 ./CarotidIso.ex
#echo Running on $(hostname --fqdn):  'Hello, world!'
