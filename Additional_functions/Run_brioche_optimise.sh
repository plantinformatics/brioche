#!/bin/bash

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="example_runbriocheoptimise"

# The project ID which this job should run under:
#SBATCH --account="default"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# The total amount of memory in megabytes in the job:
#SBATCH --mem=60GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=user123@mail.mail

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=120:0:00


module load Nextflow
module load Miniconda3
module load git

############### Brioche variables to set ######################

# this file points to the input markers file for brioche input. 
target=/filepath/briochetargetdesign.tsv

# Shortname of targets file and its function e.g., is it the targets of a SNPchip
chipname="probename"

# Location of parameter file
params="/filepath/brioche/params.config"

# absolute pathway to reference genome
genomedir="/filepathgenome/"

# shortname of reference genome 
genome="genome_1"

# fastafile of reference genoke
fasta="genome_1.fasta"

# Output directory to save brioche runs (all runs will be saved as subfolders in here)
outdir="/filepath/Brioche/testing_optimisation"

# Path to main.nf in Brioche folder (points to main.nf located in the landing folder of the brioche git 
nextflowpath="/filepath/brioche/main.nf"

# email address for updates (will give an update when each Nextflow run completes/if they fail)
emailaddress="user123@mail.mail"
####################



####### Thresholds to test ###################### 
identity_thresholds=("90" "95" "99") # pident threholds (pairwise identity to the reference genome cut off) 
coverage_thresholds=("90" "95" "97") # coverage thresholds (minimum coverage of aligned query to reference sequence cut off)
maxhits=("10") # Maximum allowed secondary hits allowed for filter secondaries stage (will remove all markers with combatable blasthits>maxhits at remove secondaries/hybrids stage [after filtering by identity and coverage])

#####################


######### Script ##################

currdir=$(pwd)

for mhit in "${maxhits[@]}"; do
  for identity in "${identity_thresholds[@]}"; do
    for coverage in "${coverage_thresholds[@]}"; do

      mkdir -p "${outdir}/${chipname}_${genome}_ident${identity}_cov${coverage}_mhits${mhit}"
      cd "${outdir}/${chipname}_${genome}_ident${identity}_cov${coverage}_mhits${mhit}"
      nextflow run "$nextflowpath" -profile slurm \
        --emailaddress "$emailaddress" \
        --genomefasta ${genomedir}/${fasta} \
        --genomename ${genome} \
        --probename ${chipname} \
        --targetdesign ${target} \
        --paramfile "$params" \
        --resultsdir "${outdir}/${chipname}_${genome}_ident${identity}_cov${coverage}_mhits${mhit}" \
        --workdir "${outdir}/${chipname}_${genome}_ident${identity}_cov${coverage}_mhits${mhit}" \
        --markercharacter "D" \
        --coverage ${coverage} \
        --pident ${identity} \
        --maximumhits ${mhit} \
        --mode prod \
        --keepduplicates false
      cd ${currdir}
    done
  done
done
