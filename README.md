# WhitefishReferenceGenome
Scripts used for Coregonus sp. "Balchen" genome assembly/analysis in De-Kayne, Zoller &amp; Feulner 2019

Parameters/submission scripts are included for:

#### Assembly using Falcon Canu and wtdbg2 
- FalconAssembly.txt - parameters for Falcon assembly by DNA Nexus
- CanuAssembly.txt - parameters for Canu assembly
- wtdbg2Assembly.txt - parameters for wtdbg2 assembly

#### Purge_haplotigs
- PurgeHaplotigs.txt - all commands for purging of haplotigs for each assembly
- submit.purge_haplotigs.sh - submission script for purge_haplotigs

#### Analysis and validation of assemblies (estimation of genome size, assessing coverage of Illumina data and confirming synteny to linkage map)
- GenomeAssessmentAnalysis.txt - all commands run for estimation of genome size and analysis of each genome assembly 
- PhaseSynteny_FR2.R - R script for coverage and synteny analysis of Falcon assembly
- PhaseSynteny_CR2.R - R script for coverage and synteny analysis of Canu assembly
- PhaseSynteny_WR2.R - R script for coverage and synteny analysis of wtdbg2 assembly

#### Repeat masking
- RepeatAnnotationCommands.txt - all commands for repeat masking of each assembly

#### Annotation
- MakerCommands.txt - all commands for MAKER2 annotation of wtdbg2 assembly

parameter files for round one of maker (remove \_r1 to use):
- maker_bopts.ctl
- maker_exe.ctl
- maker_opts_r1.ctl

modified parameter files for round 2 and round 3 of maker (remove \_r2 or \_r3 to use):
- maker_opts_r2.ctl
- maker_opts_r3.ctl








