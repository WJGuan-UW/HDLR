#!/bin/bash

### Typically, using # indicates comments in bash
# but the exact expression #SBATCH is special and not read as a comment
# the following lines indicate parameters of the slurm job that is submited

#SBATCH --partition=default_partition    # Indicates SDS partition
#SBATCH --time=4-12:00:00                # wall clock limit (d-hh:mm:ss)
#SBATCH --mem-per-cpu=5000		          # memory per cpu
#SBATCH --cpus-per-task=1		            # number of cpus per task (each array index is 1 task)
#SBATCH --mail-type=ALL			            # Send email about start/stop
#SBATCH --mail-user=wg285@cornell.edu	  # email to send to
#SBATCH --job-name=Deb_AR              # Job name
#SBATCH -o Rout/Deb_AR_%j.out         # output file (%j expands to jobID)
#SBATCH -e Rout/Deb_AR_%j.err         # error log file (%j expands to jobID)
#SBATCH --get-user-env                  # retrieve the users login environment
#SBATCH --array=1-200                   # number of jobs to submit




# This line actually submits the job
R CMD BATCH --no-save --no-restore "--args runInd=${SLURM_ARRAY_TASK_ID}" Debinf_AR.R "output/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.Rout"