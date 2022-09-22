#!/bin/bash

## specify your allocation (with the _g) and that you want GPU nodes
#SBATCH -A m4106_g
#SBATCH -C gpu

## the job will be named "Z4c2" in the queue and will save stdout to z4c2_[job ID].out
#SBATCH -J Z4c2
#SBATCH -o z4c2_%j.out

## set the max walltime
#SBATCH -t 2:00

## specify the number of nodes you want
#SBATCH -N 1

## we use the same number of MPI ranks per node as GPUs per node
#SBATCH --ntasks-per-node=4

## assign 1 MPI rank per GPU on each node
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3

# the -n argument is (--ntasks-per-node) * (-N) = (number of MPI ranks per node) * (number of nodes)
srun -n 2 ./main3d.gnu.MPI.CUDA.ex inputs
