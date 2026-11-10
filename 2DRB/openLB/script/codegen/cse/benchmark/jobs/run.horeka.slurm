#!/bin/bash
#SBATCH --partition=accelerated
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00 
#SBATCH --output=benchmark_%j.out
#SBATCH --error=benchmark_%j.err

module purge
module use /software/easybuild/modules/all/
module load OpenMPI/5.0.7-GCC-14.2.0-CUDA-12.8.0

source $HOME/roofline_pyenv/bin/activate

./run_multiple_dynamics.sh --size 25 ../target/wind.dynamics
./run_multiple_dynamics.sh --size 50 ../target/wind.dynamics
#./run_multiple_dynamics.sh --size 100 ../target/wind.dynamics
#./run_multiple_dynamics.sh --size 200 ../target/wind.dynamics
#./run_multiple_dynamics.sh --size 300 ../target/wind.dynamics


