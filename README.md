# WhitefishReferenceGenome
Scripts used for Coregonus sp. "Balchen" genome assembly/analysis in De-Kayne, Zoller &amp; Feulner 2019

For chromosome-scale scaffolds of wtdbg2 assembly see European Nucleotide Archive accession ERZ1030224 
Unscaffolded contigs (7815) are included here as two files: AWG_v1.unscaffolded1.fastq.gz and AWG_v1.unscaffolded2.fastq.gz

Parameters/submission scripts are included for:

***

### Assembly using Falcon Canu and wtdbg2 
- FalconAssembly.txt - parameters for Falcon assembly by DNA Nexus
- CanuAssembly.txt - parameters for Canu assembly
- wtdbg2Assembly.txt - parameters for wtdbg2 assembly

***

### Purge_haplotigs
- PurgeHaplotigs.txt - all commands for purging of haplotigs for each assembly
- submit.purge_haplotigs.sh - submission script for purge_haplotigs

***

### Analysis and validation of assemblies (estimation of genome size, assessing coverage of Illumina data and confirming synteny to linkage map)
- GenomeAssessmentAnalysis.txt - all commands run for estimation of genome size and analysis of each genome assembly 
  - cov.per.window.pl - for coverage analysis
  
- PhaseSynteny_FR2.R - R script for coverage and synteny analysis of Falcon assembly
  - FalconR2ChromosomeLenghts.txt - falcon assembly chromosome lengths
  - SA_linkagemap_FalconR2PhaseMapped_Filtered.csv - mapping of SA linkage map to Falcon assembly
  
- PhaseSynteny_CR2.R - R script for coverage and synteny analysis of Canu assembly
  - CanuR2ChromosomeLenghts.txt - falcon assembly chromosome lengths
  - SA_linkagemap_CanuR2PhaseMapped_Filtered.csv - mapping of SA linkage map to Canu assembly
  
- PhaseSynteny_WR2.R - R script for coverage and synteny analysis of wtdbg2 assembly
  - wtdbg2ChromosomeLenghts.txt - falcon assembly chromosome lengths
  - SA_linkagemap_Wtdbg2PhaseMapped_Filtered.csv - mapping of SA linkage map to wtdbg2 assembly
  
 - LM_stats.txt

***

### Repeat masking
- RepeatAnnotationCommands.txt - all commands for repeat masking of each assembly
- SalmonidaeRepbase.lib - Salmonidae repeat library used in repeat annotation

***

### Annotation
- MakerCommands.txt - all commands for MAKER2 annotation of wtdbg2 assembly
  - prepare.gff-merge.command.sh
  - gff_merge.command.sh
  - collect.maker.fastas.sh

parameter files for round one of maker (remove \_r1 to use):
- maker.run.cmd_r1.sh
- prep.and.start.maker.runs_r1.sh
- maker_bopts.ctl
- maker_exe.ctl
- maker_opts_r1.ctl

modified parameter files for round 2 and round 3 of maker 
(remove \_r23 to use maker.run.cmd.sh and prep.and.start.maker.runs.sh in rounds 2 and 3 of annotation)
(remove \_r2 or \_r3 to use maker_opts.ctl in rounds 2 and 3 of annotation):
- maker.run.cmd_r23.sh
- prep.and.start.maker.runs_r23.sh
- maker_opts_r2.ctl
- maker_opts_r3.ctl
