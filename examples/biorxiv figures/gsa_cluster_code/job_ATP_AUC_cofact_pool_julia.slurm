#!/bin/bash
# Job name:
#SBATCH --job-name=global_sens_analysis
#
# Account:
#SBATCH --account=fc_moorjani
#
# Partition:
#SBATCH --partition=savio
#
# Request one node:
#SBATCH --nodes=24
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=12:00:00

# Load software
module purge
module load julia/1.8.2

# Run Julia script
julia Sensitivity_analysis_ATP_AUC_cluster_cofactor_pool.jl

