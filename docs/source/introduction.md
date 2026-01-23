#Introduction

Brioche overview:
================

**Brioche** is a bioinformatics pipeline designed to help bridge the gap between static datasets and the rapidly evolving pangenomics landscape. Brioche uses Nextflow, a domain-specific language (DSL) for defining workflows, as its workflow management tool. 
Brioche streamlines marker mapping to reference genomes through a Nextflow-powered pipeline.

Brioche is funded as part of the Australian Grains Genebank Strategic Partnership, a $30M joint investment between the Victorian State Government 
and Grains Research and Development Corporation (GRDC) that aims to unlock the genetic potential of plant genetic resources for the benefit of the Australian grain growers.

In this user guide we will go over the basic structure of Brioche as well as provide in depth use cases for a variety of functionality 

Key features:
=============

**Efficient workflow management**: Nextflow, a domain-specific language designed for workflows, orchestrates processing steps with clarity and precision.
**Formalized process definition**: Each step, from input to output, is clearly defined within the pipeline, ensuring transparency and reproducibility.
**Automated job submission and connection**: Nextflow seamlessly handles job execution and data flow between steps, freeing you from manual task management.
**Datatype/species agnostic**: Brioche has been designed to work with a wide range of SNP data types including probe capture data, DArT sequencing methods, and GBS sequencing results for any species.
**Comprehensive summaries** helping users derive more clear conclusions from their data and help plan future sequencing runs.

Extensive metadata recording: Brioche effectively captures Metadata for settings used within the run and helps build VCF files containing relevant information about target reference genomes and run metadata.


Key use cases:
==============

Brioche is designed for five main usecases briefly highlighted below (see tabs for more indepth descriptions)

1. Easy remapping of markers across any distinct reference genome and the reanchoring of genotypes to the new reference genome 
{doc}`Usecase1: Remapping data across reference genomes <usecases/usecase1_remapping>`


2. The insilico genotyping of dozens + of reference genomes to generate equivalent genotype matricies for any given reference and allowing for whole genomes to be directly analysed alongside user sample genotypes
{doc}`Usecase2: In silico genotyping of reference genomes <usecases/usecase2_insilico_genotyping>`


3. Remapping of markers using a wide range of filters designed to identify which markers are hitting redundant parts of a genome and therefore are less reliable for use in certain downstream analyses
{doc}`Usecase3: Testing redundancy/accuracy in marker datasets <usecases/usecase3_redundancy_accuracy>`


4. The creation of custom marker datasets from existing analyses and testing the likely redundancy of newly designed markers against a wide range of reference genomes for a target species.
{doc}`Usecase3: Testing redundancy/accuracy in marker datasets <usecases/usecase3_redundancy_accuracy>`


5. The mapping of multiple different datasets to a unified reference genome allowing for merging across shared loci and other downstream applications (e.g., imputations)
{doc}`Usecase4: Merge datasets <usecases/usecase4_merge_datasets>`



Workflow summary:
=================

At it's core, brioche utilises BLASTn to map desired markers onto all possible positions of a given reference genome. After mapping markers, a series of layered filters are applied to identify the most probable marker position and where not determinable will position the marker as unplaced.

Brioche has five main filtering stages briefly outlined below (although, see {doc}`Setting up a run <setting_up_a_run>` for a detailed rundown of how to updated run settings and see {doc}`Citations <citations>` .


.. image:: Images/Brioche_process_only_dag.png
   :alt: Brioche dag
   :width: 300px
   :align: center


Blastn output:
~~~~~~~~~~~~~~

Markers which return no matches to a reference genome are recorded as unplaced. Several BLASTn settings can be changed from the parameters file to change the sensitivity thresholds mapping.


Identity/coverage cutoffs:
~~~~~~~~~~~~~~~~~~~~~~~~~~

Blastn hits returned from are then filtered based on the user given paiwise identity and alignment coverage % cutoffs. 


Secondary hits/hybridisation:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Blastn hits are then filtered by one-two filtered dependent on settings applied. If the user has selected settings for Probe capture design (setting istarget3primeend=TRUE) hits will be filtered to remove any which do not have the minimum set number of identical bases 
immediately before the target SNP site (setting extendablebps). Any markers with over 10 possible hits remaining at this stage are placed as not determinable


Intermediary filters:
~~~~~~~~~~~~~~~~~~~~~

The following filtering step changes depending on whether priors files are provided to inform marker position. In general, the top hits for the top chromosome, or equivalent bitscore hits across multiple chromosomes are kept.


Strict filters:
~~~~~~~~~~~~~~~

The following filtering step changes depending on whether priors files are provided to inform marker position. In general, if a marker has one clear best BLASTn hit, the marker is positioned there. If there are more than two possible positions the marker is moved to unplaced




User guide:
===========

Installation:
~~~~~~~~~~~~~

For instalation instructions go to {doc}`Installation <installation>`


Run setup:
~~~~~~~~~~

For setting up a run go to {doc}`Setting up a run <setting_up_a_run>`



Results interpretation:
~~~~~~~~~~~~~~~~~~~~~~~

For Results interpretation go to {doc}`Results <results>`



FAQ/Troubleshooting:
~~~~~~~~~~~~~~~~~~~~

For a list of FAQ/Troubleshooting questions go to {doc}`FAQ <faq>`



Citations:
~~~~~~~~~~

For a detailed list of software used go to {doc}`Citations <citations>`
