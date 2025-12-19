#!/bin/bash
set -euo pipefail


# Run example
#   bash build_brioche-vcf_env.sh

module load Miniconda3

eval "$(conda shell.bash hook)"
conda_path=$(conda info --base)
conda_base=$(conda info --base)
conda config --set channel_priority flexible

ENV_NAME="brioche-vcf"

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "Conda env '$ENV_NAME' exists proceeding with anchoring"
else
  echo "Conda env '$ENV_NAME' not found will begin download of the environment before beginning anchoring"

  mkdir -p "$conda_path"/envs/brioche-vcf

  env_path="$conda_path"/envs/brioche-vcf

  conda env create --file brioche-vcf.yaml --prefix $env_path --yes

  echo "Conda environment created at: $env_path"
  echo "Conda environment brioche-vcf fully installed"

fi
