#!/bin/bash

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="example_slurmjob"

# The project ID which this job should run under:
#SBATCH --account="default"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# The total amount of memory in megabytes in the job:
#SBATCH --mem=40GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=user123@mail.mail

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=24:0:00


#Build target table for the new brioche format
target=/filepath/briochetargetdesign.tsv

# Used in naming of results files 
chipname="probename"

# directory of reference genome
genomedir=/filepathgenome/

# Genome name without extension
genome="genomename"

# genome name with 
fasta="genomename.fa"


# param file can be set as anything but if not set in the nextflow run command below will default to the params.config present in the same directory as the main.nf
params="/filepath/brioche/params.config"
briochemain="/filepath/brioche/main.nf"

# required modules 
module load Nextflow
module load Miniconda3
module load git

    workdir="/filepath/brioche/results/${genome}_${chipname}/"
    resultsdir="/filepath/brioche/results/${genome}_${chipname}/"
    mkdir -p ${resultsdir}
    mkdir -p ${workdir}
    cd ${workdir}



    nextflow run "${briochemain}" -profile 'slurm' \
    --emailaddress user123@mail.mail \
    --genomefasta ${genomedir}/${fasta} \
    --genomename ${genome} \
    --probename ${chipname} \
    --targetdesign ${target} \
    --paramfile "${params}" \
    --resultsdir "${resultsdir}" \
    --markercharacter "D" \
    --workdir "${workdir}" \
    --coverage 70 \
    --pident 90 \
    --maximumhits 10 \
    --mode prod \
    --keepduplicates FALSE
