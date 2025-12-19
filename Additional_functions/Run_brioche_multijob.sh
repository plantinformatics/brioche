#!/bin/bash

# Multithreaded (SMP) job: must run on one node
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="example_multirun_brioche"

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
#SBATCH --time=72:0:00


module load Nextflow
module load Miniconda3
module load git

# Brioche variables to set 

# Shortname of targets file and its function e.g., is it the targets of a SNPchip
chipname="WGS"

# Location of parameter file
params="/filepath/brioche/params.config"


# absolute pathway to reference genome
genomedir="/filepathgenome/"

# shortname of reference genome (e.g., without the file extension [no .fasta, .fa, or .fasta.gz etc])
genome="genomename"

# fastafile of reference genoke
fasta="genomename.fa"

# Input directory of list of brioche target files to run Brioche on simultaneously (separate brioche instance starts up per file in this input directory)
BRIOCHE_INPUT_DIR="/filepath/brioche/brioche_inputs/"

# Output directory to save brioche runs (all runs will be saved as subfolders in here)
outdir="/filepath/brioche/Run_brioche/results"

# Path to main.nf in Brioche folder (points to main.nf located in the landing folder of the brioche git 
nextflowpath="/filepath/brioche/main.nf"

# email address for updates (will give an update when each Nextflow run completes/if they fail)
emailaddress="user123@mail.mail"

# identity cut off for filtering
identity=90
# Coverage cut off for filtering
coverage=70
# maximum number of secondary hits for filtering
mhit=10
# Partition to use 
slurm_partition="batch"
# Time per run (make sure it is longer than expected as many tasks may take longer to start in wait queues if multiple runs are being put on at the same time)
slurm_time=196:0:00
# Memory for the nextflow main script (helps with a couple of set up commands but otherwise Nextflow requests memory per task later so this doesn't need to be big
slurm_mem=40G
# CPUs for the nextflow main script (helps with a couple of set up commands but otherwise more CPUs are allocated per task later so this doesn't need to be big
slurm_cpus=4
# Set the account of the data 
slurmaccount="default"

# An initial delay is given to allow the first nextflow instance to create the software environment without having nextflow error from trying to multibuild
INITIAL_SUBMIT_DELAY=900

####################

currdir=$(pwd)
mkdir -p "${outdir}"

# Collect inputs (adjust the pattern if needed)
mapfile -t INPUTS < <(find "${BRIOCHE_INPUT_DIR}" -maxdepth 1 -type f -name 'Brioche_inputfile_*.tsv' | sort)

if (( ${#INPUTS[@]} == 0 )); then
  echo "No Brioche input files found in: ${BRIOCHE_INPUT_DIR}"
  exit 1
fi

echo "Submitting ${#INPUTS[@]} Brioche runs…"

i=0
for target in "${INPUTS[@]}"; do
  bname="$(basename "$target" .tsv)"
  chrom="${bname#Brioche_inputfile_}"
  run_dir="${outdir}/${chipname}_${genome}_${chrom}_ident${identity}_cov${coverage}_mhits${mhit}"
  mkdir -p "${run_dir}"

  jobid=$(
    sbatch --parsable <<SBATCH_EOF
#!/bin/bash
#SBATCH -J brioche_${chrom}
#SBATCH --account=${slurmaccount}
#SBATCH -p ${slurm_partition}
#SBATCH -t ${slurm_time}
#SBATCH -c ${slurm_cpus}
#SBATCH --mem=${slurm_mem}
#SBATCH -o ${run_dir}/slurm-%j.out
#SBATCH -e ${run_dir}/slurm-%j.err
set -euo pipefail
mkdir -p "${run_dir}"
cd "${run_dir}"
module load Nextflow
module load Miniconda3
module load git
nextflow run "${nextflowpath}" -profile slurm -resume \
  --emailaddress "${emailaddress}" \
  --genomefasta "${genomedir}/${fasta}" \
  --genomename "${genome}" \
  --probename "${chipname}" \
  --targetdesign "${target}" \
  --paramfile "${params}" \
  --resultsdir "${run_dir}" \
  --workdir "${run_dir}" \
  --markercharacter "D" \
  --coverage "${coverage}" \
  --pident "${identity}" \
  --maximumhits "${mhit}" \
  --mode prod \
  --keepduplicates false
SBATCH_EOF
  )

  submit_time=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$submit_time] Submitted ${bname} (chrom=${chrom}) ? job ${jobid}"

  if (( i == 0 )) && (( ${#INPUTS[@]} > 1 )) && (( INITIAL_SUBMIT_DELAY > 0 )); then
    echo "Waiting ${INITIAL_SUBMIT_DELAY}s before submitting remaining $(( ${#INPUTS[@]} - 1 )) jobs…"
    sleep "${INITIAL_SUBMIT_DELAY}"
  fi
  ((i++))
done

echo "All jobs submitted."
cd "${currdir}"
