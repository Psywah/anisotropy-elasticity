#SBATCH --job-name mpi_test
#SBATCH --qos janus
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 12
#SBATCH --time 00:10:00
#SBATCH --output mpi_test.out

# the slurm module provides the srun command
module load slurm

module load intel/impi-13.0.0

srun mpi_test
