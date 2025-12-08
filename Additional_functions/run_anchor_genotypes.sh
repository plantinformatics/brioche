#!/bin/bash

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="example_anchoring"

# The project ID which this job should run under:
#SBATCH --account="default"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

# The total amount of memory in megabytes in the job:
#SBATCH --mem=60GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=user123@mail.mail

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=4:0:00

module load R/

date

Rscript Anchoring_script.R --Rawgenotypes Rawgenotypesfile.tsv --Briochemappings Briocheresultsdatasetx_all_markers1to1stagingforvcf.csv --IsVCFraws false --IsSNPchip true --isDArTfile false --Outputfilename "dataset1_anchored_genotypes.vcf" --genomename "Genome123XYZ" --Accession "GCF000XYZ"

date



