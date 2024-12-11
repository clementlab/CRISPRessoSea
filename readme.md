# CRISPRessoSea
CRISPRessoSea is a tool for processing genome editing from multiple amplicon sequencing experiments.

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
```
usage: CRISPRessoSea.py Process [-h] [-o OUTPUT_FOLDER] -g GUIDE_FILE -s SAMPLE_FILE -x GENOME_FILE [-p N_PROCESSES] [--skip_bad_chrs] [--plot_only_complete_guides] [--min_amplicon_coverage MIN_AMPLICON_COVERAGE]
                                [--sort_based_on_mismatch] [--allow_guide_match_to_other_region_loc]

options:
  -h, --help            show this help message and exit
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder
  -g GUIDE_FILE, --guide_file GUIDE_FILE
                        Guide file - list of guides with one guide per line
  -s SAMPLE_FILE, --sample_file SAMPLE_FILE
                        Sample file - list of samples with one sample per line
  -x GENOME_FILE, --genome_file GENOME_FILE
                        Bowtie2-indexed genome file - files ending in .bt2 must be present in the folder.
  -p N_PROCESSES, --n_processes N_PROCESSES
                        Number of processes to use
  --skip_bad_chrs       Skip regions on bad chromosomes (chrUn, random, etc)
  --plot_only_complete_guides
                        Plot only guides with all values. If not set, all guides will be plotted.
  --min_amplicon_coverage MIN_AMPLICON_COVERAGE
                        Minimum number of reads to cover a location for it to be plotted. Otherwise, it will be set as NA
  --sort_based_on_mismatch
                        Sort guides based on mismatch count. If true, the on-target will always be first
  --allow_guide_match_to_other_region_loc
                        If true, guides can match to regions even if the guide chr:start is not in that region (e.g. if the guide sequence is found in that region). If false/unset, guides can only match to regions matching the guide chr:start position. This flag should be set if the genome for guide design was not the same as the analysis genome.

```

### Replot
```
usage: CRISPRessoSea.py Replot [-h] [-o OUTPUT_FOLDER] [-p FILE_PREFIX] -f REORDERED_GUIDE_FILE [-n NAME_COLUMN]

options:
  -h, --help            show this help message and exit
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder
  -p FILE_PREFIX, --file_prefix FILE_PREFIX
                        File prefix for output files
  -f REORDERED_GUIDE_FILE, --reordered_guide_file REORDERED_GUIDE_FILE
                        Reordered guide file - made by reordering rows from aggregated_stats_all.txt
  -n NAME_COLUMN, --name_column NAME_COLUMN
                        Column name to set as the displayed name for each sample in the plot
```


