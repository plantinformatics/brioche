# Brioche
Brioche is a bioinformatics pipeline for mapping markers to a given reference genome. brioche uses [Nextflow](https://www.nextflow.io/), a domain-specific language (DSL) for defining
workflows, as its workflow management tool. Brioche streamlines marker mapping to reference genomes through a Nextflow-powered pipeline.

Brioche is funded as part of the Australian Grains Genebank Strategic Partnership, a $30M joint investment between the Victorian State Government and Grains Research and Development Corporation (GRDC) that aims to unlock the genetic potential of plant genetic resources for the benefit of the Australian grain growers.


**Key features:**
- Efficient workflow management: Nextflow, a domain-specific language designed for workflows, orchestrates processing steps with clarity and precision.
- Formalized process definition: Each step, from input to output, is clearly defined within the pipeline, ensuring transparency and reproducibility.
- Automated job submission and connection: Nextflow seamlessly handles job execution and data flow between steps, freeing you from manual task management.

## Getting started
To get started, you will need to load or install **Nextflow**.

* Load nextflow using the 'module load' option if running on basc
  ```{bash}
  module load Nextflow 
  ```

* Alternatively, you can install Nextflow by using the following command (optional): 
  ```{bash}
  wget -qO- https://get.nextflow.io | bash 
  ```
  see [Getting started with nextflow](https://www.nextflow.io/docs/latest/getstarted.html) for further details.
  
* Download brioche and then cd to the brioche directory : 

  ```{bash}
  git clone https://github.com/plantinformatics/brioche.git
  cd brioche
  ```

## Running brioche with default params

* Launch the pipeline execution with the following command:

**local** 
  ```{bash}
  nextflow run main.nf --mode "test"
  ```
**Using slurm**
  ```
  nextflow run main.nf -profile 'slurm' --mode "test"
  ```

Note that this runs the pipeline with the default parameters and test data in [*Data*](Data);

## Running brioche with new params

To run with new parameters, edit the [*params.config*](params.config) file and then run ;
  ```{bash} 
  nextflow run main.nf --mode 'prod' --paramfile 'path 2 params.config file'
  ```
alternatively, parameters can also be passed via the commandline

  ```{bash}
  nextflow run main.nf [options]
  Options:
  --mode set to 'prod'
  --genomefasta Absolute path to reference genome fasta file
  --genomename Name of reference genome
  --probename Name of probe
  --targetdesign Absolute path to target file of the probe
  --markercharacter Character used to replace target marker in probe sequence
  ```
Note the format of the target design table should be a tab delimited file with the following columns;


| Column | Description |
| -- | -- |
| 1. ID  | The unique ID for each probe  |
| 2. Probe Sequence  | The sequence of each probe |
| 3. Target.bp   | Target position of bp in the probe sequence |
| 4. Target.base  | Details about the Targetted base |


See [*example target design table*](Data/AVRGRDC_Pulses_v1_20006795X370754_A2_Chickpea-target-new-format.tsv)



