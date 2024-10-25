#!/bin/bash

#SBATCH --job-name=simulation-scenario3  # Job name
#SBATCH --array=1-100     # Number of jobs
#SBATCH --partition=standard    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/simulation-scenario3/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/simulation-scenario3/slurm-%x-%a-%J.out # Same file for output and error

R CMD BATCH --vanilla BATCH_simulation-scenario3.R output/simulation-scenario3/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd c:/Users/hpl802/Documents/Github/multiverse_tools/sensitivity_analysis/
                      ## cd /projects/biostat01/people/hpl802/pipeline/
## clean outputs      ## rm -r ./output/simulation-scenario3/*
## clean results      ## rm -r ./Results/simulation-scenario3/*
## submission command ## sbatch SUBM_simulation-scenario3.slurm

## submission output  ## Submitted batch job 2682252 (time using C-c t: torsdag 22-12-08 at 12:10) 

## init               ## module module load gcc/11.2.0 R/4.1.2
                      ## cat .bash_profile (in the root directory on the server to visualize the current configuration)
                      ## vi .bash_profile (in the root directory on the server to modify the current configuration)

## commands           ## squeue           : view current jobs 
                      ## squeue -u id     : view current jobs for user id
                      ## scancel          : delete the job with job id 1034 from the queue type
                      ## sinfo            : view state of the queues
                      ## sshare           : check own usage of the server
                        
## documentation      ## from SLURM: https://slurm.schedmd.com/pdfs/summary.pdf: key functions
                      ##           : https://slurm.schedmd.com/sbatch.html (filename pattern)
                      ## from KU   : https://hpc.ku.dk/documentation/slurm.html
                      ##           : https://github.com/aejensen/computing-biostat-ucph/blob/main/README.md
