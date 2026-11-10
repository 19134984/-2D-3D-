#!/bin/bash 
#SBATCH --partition=dev-g 
#SBATCH --nodes=1 
#SBATCH --exclusive
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=7 
#SBATCH --gpus-per-node=1
#SBATCH --time=01:00:00 
#SBATCH --account=project_XXXX
#SBATCH --output=benchmark_%j.out 
#SBATCH --error=benchmark_%j.err 
#SBATCH --mem=0

module --force purge > /dev/null 2>&1 
module load LUMI/25.03 partition/G PrgEnv-amd rocm cray-python

SIZE=$(./babel_size.sh $1)

cd BabelStream

echo "Size: ${SIZE}"

./build/hip-stream --float -s $SIZE

SCRIPT_HOME=$OLB_ROOT/script/codegen/cse

cd ..

# PYTHON ENVIRONMENT MUST BE INSTALLED IN ADVANCE USING THE $SCRIPT_HOME/pyenv_requirements.txt 
source ../.roofline_pyenv/bin/activate

./run_multiple_dynamics.sh --size $1 ../target/roofline.dynamics
