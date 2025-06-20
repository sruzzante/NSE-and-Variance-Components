#!/bin/bash -l
#SBATCH --job-name run_nh
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-task=2
#SBATCH --mem=15G 
##SBATCH --gres=gpu:1
#SBATCH --mail-user=sruzzante@uvic.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load python
module load netcdf
module load scipy-stack/2023b


source ../../NH_ENV4/bin/activate

#python run_nh-ens0.py
#python run_nh-ens1.py
python run_nh-ens2.py
#python run_nh-ens3.py
#python run_nh-ens4.py
