# CRISPRessoOcean
CRISPRessoOcean is a tool for processing genome editing from multiple amplicon sequencing experiments.

The main program has two subcommands:
### Process
`Process`: This program will process a pool of pools given:
 -  a guide information file with headers: Name (Name of the on-target sequence), Sequence, PAM (PAM sequence), Score (optional), #MM (Number of mismatches), Gene (optional), Locus (formatted like chr1:+2345), and anno (optional, if provided, this name will be shown)
 -  a sample file with headers: Name, fastq_r1, fastq_r2 (optional for single-end reads)
 -  a reference genome. The path to the reference genome and bowtie2 indices is provided not including the trailing .fa

 ### Replot
 `Replot`: This program will replot data from a finished analysis performed using `Process`. This avoids reprocessing or rerunning any analyses but allows users to reorder or subset guides to be plotted. This program accepts as input:
 - A modified file derived from an 'aggregated_stats_all.txt' output from a completed `Process` run.


## Complete command description:
### Process
Example: 



