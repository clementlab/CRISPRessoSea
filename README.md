# CRISPRessoSea
CRISPRessoSea is a tool for processing genome editing from multiple pooled amplicon sequencing experiments to characterizing editing at on- and off-targets across experimental conditions. The tool accepts raw sequencing files as input, performs analysis at each site in each sample, performs statistical analysis of editing, and produces plots and reports detailing editing rates. 

## Installation
CRISPRessoSea can be installed into an environment that contains CRISPResso2 and its dependencies.
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -y -n CRISPRessoSea crispresso2 # Create CRISPRessoSea conda environment
conda activate CRISPRessoSea
pip install git+https://github.com/clementlab/CRISPRessoSea.git
```

## Running CRISPRessoSea

CRISPRessoSea operates in three primary running modes to streamline the preparation, processing, and replotting of large-scale pooled sequencing experiments to measure and compare CRISPR genome editing.

### MakeGuideFile
`MakeGuideFile`: This mode enables users to generate a comprehensive target file that includes both on- and off-target sites for one or more guide RNA sequences. Off-target sites are computationally predicted using [Cas-OFFinder](https://github.com/snugel/cas-offinder), which reports genomic coordinates, sequences, and mismatch counts for each potential target. This functionality is especially useful for designing pooled experiments to profile editing activity at predicted off-target sites, or for preparing analyses when the off-target sequences or locations are not already known. The output of this mode is a standardized target file compatible with CRISPRessoSea’s Process mode.

### Process
`Process`: This program will process a pool of pools given:
 -  a target information file with headers: Guide (Name of the on-target sequence guide), Sequence, PAM (PAM sequence), #MM (Number of mismatches), and Locus (formatted like chr1:+2345). Target names can be provided in a 'Target' column. 
 -  a sample file with headers: Name, fastq_r1, fastq_r2 (optional for single-end reads). A column specifying the group can also be provided in this file.
 -  a reference genome. The path to the reference genome and bowtie2 indices is provided not including the trailing .fa

 ### Replot
 `Replot`: This program will replot data from a finished analysis performed using `Process`. This avoids reprocessing or rerunning any analyses but allows users to reorder or subset guides to be plotted. This program accepts as input:
 - A modified file derived from an 'aggregated_stats_all.txt' output from a completed `Process` run.


## Complete command description:
### MakeGuideFile
```
usage: CRISPRessoSea.py MakeGuideFile [-h] [-o OUTPUT_FOLDER] [-p FILE_PREFIX] -g GUIDE_SEQ [-gn GUIDE_NAME] [--pam PAM] -x GENOME_FILE [--max_mismatches MAX_MISMATCHES] [--max_dna_bulges MAX_DNA_BULGES] [--max_rna_bulges MAX_RNA_BULGES] [-v VERBOSITY] [--debug]

options:
  -h, --help            show this help message and exit
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder (default: None)
  -p FILE_PREFIX, --file_prefix FILE_PREFIX
                        File prefix for output files (default: CRISPRessoSea)
  -g GUIDE_SEQ, --guide_seq GUIDE_SEQ
                        Guide sequence(s) to create a guide file for. Multiple guides may be separated by commas. (default: None)
  -gn GUIDE_NAME, --guide_name GUIDE_NAME
                        Guide name(s) to use for the guide sequence(s). Multiple names may be separated by commas. (default: guide_0)
  --pam PAM             PAM sequence(s) to use for the guide sequence(s). (default: NGG)
  -x GENOME_FILE, --genome_file GENOME_FILE
                        Bowtie2-indexed genome file - files ending in and .bt2 must be present in the same folder. (default: None)
  --max_mismatches MAX_MISMATCHES
                        Maximum number of mismatches to allow in the discovered offtargets (default: 4)
  --max_dna_bulges MAX_DNA_BULGES
                        Maximum number of DNA bulges to allow in the discovered offtargets. Note that Cas-OFFinder 3 is required to detect sites with bulges. (default: 0)
  --max_rna_bulges MAX_RNA_BULGES
                        Maximum number of RNA bulges to allow in the discovered offtargets. Note that Cas-OFFinder 3 is required to detect sites with bulges. (default: 0)
  -v VERBOSITY, --verbosity VERBOSITY
                        Verbosity level of output to the console (1-4) 4 is the most verbose (default: 3)
  --debug               Print debug information (default: False)
```

### Process
```
usage: CRISPRessoSea.py Process [-h] [-o OUTPUT_FOLDER] -t TARGET_FILE -s SAMPLE_FILE [--gene_annotations GENE_ANNOTATIONS] -x GENOME_FILE [-p N_PROCESSES] [--crispresso_quantification_window_center CRISPRESSO_QUANTIFICATION_WINDOW_CENTER]
                                [--crispresso_quantification_window_size CRISPRESSO_QUANTIFICATION_WINDOW_SIZE] [--crispresso_base_editor_output] [--crispresso_default_min_aln_score CRISPRESSO_DEFAULT_MIN_ALN_SCORE] [--crispresso_plot_window_size CRISPRESSO_PLOT_WINDOW_SIZE]
                                [--allow_unplaced_chrs] [--plot_only_complete_targets] [--min_amplicon_coverage MIN_AMPLICON_COVERAGE] [--sort_based_on_mismatch] [--allow_target_match_to_other_region_loc] [--top_percent_cutoff TOP_PERCENT_CUTOFF]
                                [--min_amplicon_len MIN_AMPLICON_LEN] [--fail_on_pooled_fail] [--plot_group_order PLOT_GROUP_ORDER] [--sig_method_parameters SIG_METHOD_PARAMETERS] [-v VERBOSITY] [--debug]

options:
  -h, --help            show this help message and exit
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder (default: None)
  -t TARGET_FILE, --target_file TARGET_FILE
                        Target file - list of targets with one target per line and the headers ['Guide','Sequence','PAM','#MM','Locus'] (default: None)
  -s SAMPLE_FILE, --sample_file SAMPLE_FILE
                        Sample file - list of samples with one sample per line with headers ['Name','fastq_r1','fastq_r2'] (default: None)
  --gene_annotations GENE_ANNOTATIONS
                        Gene annotations file - a tab-separated .bed file with gene annotations with columns "chr" or "chrom", "start" or "txstart", and "end" or "txend" as well as "name" (default: None)
  -x GENOME_FILE, --genome_file GENOME_FILE
                        Bowtie2-indexed genome file - files ending in and .bt2 must be present in the same folder. (default: None)
  -p N_PROCESSES, --n_processes N_PROCESSES
                        Number of processes to use. Set to "max" to use all available processors. (default: 8)
  --crispresso_quantification_window_center CRISPRESSO_QUANTIFICATION_WINDOW_CENTER
                        Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is
                        suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17 to only include mutations near the 5' end of the sgRNA.
                        (default: -3)
  --crispresso_quantification_window_size CRISPRESSO_QUANTIFICATION_WINDOW_SIZE
                        Size (in bp) of the quantification window extending from the position specified by the '--cleavage_offset' or '--quantification_window_center' parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the
                        quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp.
                        (default: 1)
  --crispresso_base_editor_output
                        Outputs plots and tables to aid in analysis of base editor studies. (default: False)
  --crispresso_default_min_aln_score CRISPRESSO_DEFAULT_MIN_ALN_SCORE
                        Default minimum homology score for a read to align to a reference amplicon. (default: 60)
  --crispresso_plot_window_size CRISPRESSO_PLOT_WINDOW_SIZE
                        Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the quantification_window_center for each guide are plotted. (default: 20)
  --allow_unplaced_chrs
                        Allow regions on unplaced chromosomes (chrUn, random, etc). By default, regions on these chromosomes are excluded. If set, regions on these chromosomes will be included. (default: False)
  --plot_only_complete_targets
                        Plot only targets with all values. If not set, all targets will be plotted. (default: False)
  --min_amplicon_coverage MIN_AMPLICON_COVERAGE
                        Minimum number of reads to cover a location for it to be plotted. Otherwise, it will be set as NA (default: 10)
  --sort_based_on_mismatch
                        Sort targets based on mismatch count. If true, the on-target will always be first (default: False)
  --allow_target_match_to_other_region_loc
                        If true, targets can match to regions even if the target chr:start is not in that region (e.g. if the target sequence is found in that region). If false/unset, targets can only match to regions matching the target chr:start position. This flag should be set
                        if the genome for guide design was not the same as the analysis genome. (default: False)
  --top_percent_cutoff TOP_PERCENT_CUTOFF
                        The top percent of aligned regions (by region read depth) to consider in finding non-overlapping regions during demultiplexing. This is a float between 0 and 1. For example, if set to 0.2, the top 20% of regions (by read depth) will be considered. (default:
                        0.2)
  --min_amplicon_len MIN_AMPLICON_LEN
                        The minimum length of an amplicon to consider in finding non-overlapping regions during demultiplexing. Amplicons shorter than this will be ignored. (default: 50)
  --fail_on_pooled_fail
                        If true, fail if any pooled CRISPResso run fails. By default, processing will continue even if sub-CRISPResso commands fail. (default: False)
  --plot_group_order PLOT_GROUP_ORDER
                        Order of the groups to plot (if None, the groups will be sorted alphabetically) (default: None)
  --sig_method_parameters SIG_METHOD_PARAMETERS
                        Parameters for the significance method in the form of: none hard_cutoff,cutoff neg_binomial,group1,group2,alpha t_test,group1,group2,alpha mean_diff,group1,group2,cutoff (default: None)
  -v VERBOSITY, --verbosity VERBOSITY
                        Verbosity level of output to the console (1-4) 4 is the most verbose (default: 3)
  --debug               Print debug information (default: False)
```

### Replot
```
usage: CRISPRessoSea.py Replot [-h] [-o OUTPUT_FOLDER] [-p FILE_PREFIX] -f REORDERED_STATS_FILE -s REORDERED_SAMPLE_FILE [-n NAME_COLUMN] [--fig_width FIG_WIDTH] [--fig_height FIG_HEIGHT] [--seq_plot_ratio SEQ_PLOT_RATIO] [--plot_group_order PLOT_GROUP_ORDER]
                               [--sig_method_parameters SIG_METHOD_PARAMETERS] [--dot_plot_ylims DOT_PLOT_YLIMS] [--heatmap_max_value HEATMAP_MAX_VALUE] [--heatmap_min_value HEATMAP_MIN_VALUE] [-v VERBOSITY] [--debug]

options:
  -h, --help            show this help message and exit
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder where output plots should be placed (default: None)
  -p FILE_PREFIX, --file_prefix FILE_PREFIX
                        File prefix for output files (default: CRISPRessoSea)
  -f REORDERED_STATS_FILE, --reordered_stats_file REORDERED_STATS_FILE
                        Reordered statistics file - made by reordering rows from aggregated_stats_all.txt (default: None)
  -s REORDERED_SAMPLE_FILE, --reordered_sample_file REORDERED_SAMPLE_FILE
                        Reordered_sample_file - path to the sample file with headers: Name, group, fastq_r1, fastq_r2 (group is always optional, fastq_r2 is optional for single-end reads) (default: None)
  -n NAME_COLUMN, --name_column NAME_COLUMN
                        Column name to set as the displayed name for each sample in the plot (default: None)
  --fig_width FIG_WIDTH
                        Width of the figure (default: 24)
  --fig_height FIG_HEIGHT
                        Height of the figure (default: 24)
  --seq_plot_ratio SEQ_PLOT_RATIO
                        Ratio of the width of the sequence plot to the data plot (>1 means the seq plot is larger than the data plot) (default: 1)
  --plot_group_order PLOT_GROUP_ORDER
                        Order of the groups to plot (if None, the groups will be sorted alphabetically) (default: None)
  --sig_method_parameters SIG_METHOD_PARAMETERS
                        Parameters for the significance method in the form of: none hard_cutoff,cutoff neg_binomial,group1,group2,alpha t_test,group1,group2,alpha mean_diff,group1,group2,cutoff (default: None)
  --dot_plot_ylims DOT_PLOT_YLIMS
                        Comma-separated min,max y-axis limits for the dot plot. If None, the y-axis limits will be set automatically. (default: None,None)
  --heatmap_max_value HEATMAP_MAX_VALUE
                        Maximum value for the heatmap color scale, where a value of 1 sets the max value color to 1% (if None, the maximum value will be determined automatically) (default: None)
  --heatmap_min_value HEATMAP_MIN_VALUE
                        Minimum value for the heatmap color scale, where a value of 1 sets the min value color to 1% (if None, the minimum value will be determined automatically) (default: None)
  -v VERBOSITY, --verbosity VERBOSITY
                        Verbosity level of output to the console (1-4) 4 is the most verbose (default: 3)
  --debug               Print debug information (default: False)
```
