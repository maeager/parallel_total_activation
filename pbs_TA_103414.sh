#!/bin/bash
# Set max wallclock time to one minute
#PBS -l walltime=4:00:00:00

#PBS -A NCIdw3

# Set name of job shown in showq
#PBS -N pbs_TotalActivation_103414

# Use submission environment
#PBS -V

# Run the job on one node, with one processor per node and one GPU per node
#PBS -l nodes=1:ppn=10
#PBS -l pmem=4000MB

# Start job from the directory it was submitted
cd /gpfs/M2Home/eagerm/Monash016/TA/singleparTotalActivation #$PBS_O_WORKDIR

#If you need to load the environment
module load matlab/r2011a

# Run a commmand that shows where the job ran
hostname

# Lets also look at what PBS ennvironment variables our environment has
env | grep PBS

time matlab -nosplash -nodesktop -r "matlabpool close force local; matlabpool local 8; fastActivation('../HCP_10_Data/103414'); matlabpool close force; exit"

echo "Job completed ok!"
