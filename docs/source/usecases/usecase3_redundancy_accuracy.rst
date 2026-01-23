Usecase3: Testing redundancy/accuracy in marker datasets
========================================================

This usecase describes how to use Brioche to test the redundancy/ accuracy of either existing or while designing marker datasets


Requirements 
~~~~~~~~~~~~

1. Brioche targets file 

2. Reference genome FASTA/FASTAs

3. An updated params.config file ready to run Brioche


Runfile
~~~~~~~

To test for the redundancy/accuracy of markers in different settings a simple loop script is provided at Additional_functions/Run_brioche_optimise.sh
This file will allow the user to set multiple values for four variables to test as combinations to investigate how the proportion of uniquely mapping markers changes

Additional instructions are provided in the script but the main components to change are

1) The Slurm ``sbatch`` header commands (e.g., ``#SBATCH --job-name="example_runbriocheoptimise"``)

2) The different variable settings:

   a) Brioche-specific settings (lines 34-63, 9 settings).

   b) Update the method for calling the required software (Nextflow, conda, and git) (lines 64-71). These are currently written as ``module load`` statements.

   c) Thresholds to test (lines 73-79, 4 variable types).



Output
~~~~~~

Brioche will be run for each combination of thresholds e.g., if identity_thresholds coverage_thresholds, maxhits, and wordsize 
all have 3 values each 81 (3^4) runs of brioche will be done so don't add too many at once to test.   

Brioche outputs will be placed in the chosen output folder with standard Brioche results but with
the name of each folder inlcuding the variables run in it e.g., Brioche_run_ident90cov70wordsize13


Extensions
~~~~~~~~~~

If the user is interested in identifying the most strictly mapped markers to design future markers from, this script allows for exploring what are the most strict parameters for a given marker dataset

After identifying strict settings, the user can then anchor the results as described in usecase 1 :doc:`Usecase1: Remapping data across reference genomes <usecase1_remapping>`
take the VCF of only mapped reads and extract markers from it using the Additional_functions/run_convert_vcf_to_brioche_input.sh
and then run Brioche on a new reference genome iteratively restricting the dataset to only the best mapping markers across all reference genomes



Other usecases
--------------

If you are interested in other usecases see.

1. If you are interested in easy remapping of markers across any distinct reference genome and the reanchoring of genotypes to the new reference genome 
:doc:`Usecase1: Remapping data across reference genomes <usecase1_remapping>`

2. If you are interested in extracting genotype calls from one or multiple reference genomes and adding them to your population genomics study go to 
:doc:`Usecase2: In silico genotyping of reference genomes <usecase2_insilico_genotyping>`

5. If you are interested in the mapping of multiple different datasets to a unified reference genome allowing for merging across shared loci and other downstream applications (e.g., imputations)
:doc:`Usecase4: Merge datasets <usecase4_merge_datasets>`


otherwise, to return to the Introduction page go to :doc:`Introduction <../introduction>`


