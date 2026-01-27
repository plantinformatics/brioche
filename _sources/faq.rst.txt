FAQ
===

1. How many markers can brioche analyse in a single analysis?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Brioche is capable of handling large numbers of markers through splitting some of the more compute and memory intensive steps in the pipeline up, but some steps require the holding of complete datasets in memory or reading/writing large files.
These steps can greatly slow down Brioche when too many markers are analysed and greatly increase the required memory to complete. From testing, Brioche starts to significantly slow down once 1 million markers are analysed at once and we recommend 
not running anything more than 2 million markers at a time. For larger datasets (large whole genome sequencing) we recommend that the markers are split by chromosome in the input file e.g., if there are 20 chromosomes and 40 million markers, splitting by
chromosome will lead to an average of ~2 million markers per chromosome and per run. Brioche can then be run on each chromosome separately, and then the 1:1 mappings files can be concatenated together for downstream applications. 
Other ways to greatly speed up Brioche when working with large numbers of markers include setting stricter blastn settings by passing a -pident and -coverage filter to blastn using the parameter otherblastoptions in the params.config file
For users who have the compute resources to run multiple Brioche runs concurrently a prewritten multijob launch script is avaliable at Additional_functions/Run_brioche_multijob.sh. This can be combined with splitting datasets by chromosome for maximum efficiency and speed


2. What are the best settings to run Brioche on?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best settings to run Brioche on will depend on both the goal of the project and the type of sequence data and so there are no best settings to run Brioche on. Generally though, if you are trying to run analyses which require exact marker positions e.g., imputation stricter is better. 
If you are running datasets generated from data like probe capture, try to set settings which reflect the probe sensitivity e.g., if the probe fails to bind at <95% sequence similarity that should be the identity threshold.
If you are running datasets generated from Reduced representation sequencing lower identity thresholds should be considered to reflect the broad capability of the sequencing methods to sequence broadly similar sequence cut sites. 


3. The DArT scripts don't seem to work?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DArT scripts to either create the Brioche input targets file or anchor the genotypes to new genomes are designed to only work with 1-row format DArT reports. This means that if the report you have is in 2-row format the scripts won't work. 
If your format is 1-row and the scripts aren't working, let us know, there is a lot of variability in what is provided in a DArT report format so we might need to update the scripts to work against different formats we didn't see during testing.


4. Some analysis steps seem to be taking a very long time or are restarting why is that?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If some steps are taking a really long time to run or a restarting it might be an indication that not enough memory is allocated to the steps. The memory settings for different steps can be changed in the file nextflow.config in the main Brioche directory.
Try increasing the memory for the section which is withLabel: 'large' or the section withLabel 'merge'. If this doesn't help or you do not have access to more compute resources, you might need to run brioche in subsets e.g., by splitting your markers by chromosome


5. Brioche seems to crash and not install the conda environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some HPC systems require the user to download/install programs in a special download/install partition of the HPC. If this is the case, Brioche may fail is just run normally. We recommend that you launch brioche (e.g., the example_run_script_brioche.sh script) with your install/download partition selected.
For the additional functions which require the conda environment brioche-vcf being built (e.g., anchoring genotypes, or converting DArT/VCF formats to brioche input) install the conda environment first by running the script build_brioche-vcf_env.sh in the Additional_functions folder on the download/install partition of your HPC


6. I don't have access to many modules and can't load Nextflow/R etc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't have access to any of the preinstalled modules we recommend you begin by trying to install Miniconda following the instructions outlined in the installation section of this manual. After that, you can try and install R, Nextflow and Git into your 
conda base environment using conda install. From there, you will need to add in a few extra lines to each brioche launch script to make sure conda is being called before you begin e.g.,

   .. code-block:: bash

      eval "$(conda shell.bash hook)"

      conda activate base


7. My HPC does not use slurm it uses PBS/QSUB/LSF/SGE etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently Brioche is written to work in a slurm management system but it may be possible for an experienced user to change the required settings to match the requirements of different HPC management systems. 
We have not tested this but feel free to try and update the ##SBATCH commands and the settings in the nextflow.config file to reflect your HPC management system. Feel free to let us know how it works out. If there is enough demand we 
might look into ways to make Brioche more HPC management system agnostic.


8. Brioche crashed and nothing above helps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are still actively testing and building on Brioche. If you experience any unusual behaviour or crashes feel free to reach out or create an issue in the github. We appreciate the feedback and any fixed bug is a step towards Brioche working perfectly.



