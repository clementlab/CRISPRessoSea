# CRISPRessoSea
CRISPRessoSea is a tool for processing genome editing from multiple pooled amplicon sequencing experiments to characterizing editing at on- and off-targets across experimental conditions. The tool accepts raw sequencing files as input, performs analysis at each site in each sample, performs statistical analysis of editing, and produces plots and reports detailing editing rates. 

[Installation](#installation) | [Running modes](#running-crispressosea) | [Tutorial](#tutorial) | [Parameters](#complete-command-description)

## Installation
CRISPRessoSea can be installed into an environment that contains CRISPResso2 and its dependencies. cas-offinder is optional but required for creating a guide info file from scratch.
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -y -n CRISPRessoSea bioconda::crispresso2 bioconda::cas-offinder # Create CRISPRessoSea conda environment
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

## Tutorial

### Download tutorial dataset
Download the tutorial dataset and change into that directory by running:
```
wget https://github.com/clementlab/CRISPRessoSea/raw/refs/heads/main/demo/small_demo.tar.gz
tar -xzf small_demo.tar.gz
cd CRISPRessoSea_demo
```
This tutorial includes a subset of data from [Cicera et al. 2020](https://www.nature.com/articles/s41587-020-0555-7) investigating the CTLA4_Site9 guide with sequence GGACTGAGGGCCATGGACACNGG. The complete dataset can be found at [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA625995](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA625995).

For ease of use, I've [created a super-small genome](demo/make_demo/small_demo/01_make_demo.ipynb) that includes the genomic sequence of only the on-target and three off-by-2 off-targets. The genomic sequence is contained in 'demo_genome.fa' and bowtie2 indices for the genome were generated using `bowtie2-build demo_genome.fa demo_genome`. 

### MakeGuideFile
 You can create a guide info file (e.g. for designing a pooled experiment to profile off-targets) by running:

```
CRISPRessoSea MakeGuideFile --guide_seq GGACTGAGGGCCATGGACAC --pam NGG --guide_name CTLA4_site9 --max_mismatches 2 --genome_file demo_genome.fa
```
- The `--guide_seq` parameter specify the on-target guide sequence and does not include the PAM sequence.
- The `--guide_name` is for convenience, and all off-targets are annotated with this guide name (see `Guide` column in output below).
- The `--pam` and `--max_mismatches` parameters are used for finding off-target locations. Here, for our small example we'll search for up to 2 mismatches, but in practice up to 4 or 5 mismatches are investigated.
- The `--genome_file` parameter specifies the path to the genome file. Here we are using a super-small genome with only the on-target and three off-by-2 off-targets.
Note that this requres [Cas-offinder](https://github.com/snugel/cas-offinder) for enumerating off-target sites.

This produces a guide info file `CRISPRessoSea_MakeGuideFileOutput//CRISPRessoSea.guide_info.txt` that contains the on- and off-target locations for the CTLA4_site9 guide.:
```
Guide   Target  Sequence        PAM     #MM     Locus   Mismatch_info
CTLA4_site9     CTLA4_site9_ON_CTLA4_site0_500  GGACTGAGGGCCATGGACAC    GGG     0       CTLA4_site0:+500        Mismatches: 0
CTLA4_site9     CTLA4_site9_OB2_CTLA4_site1_500 GGACaGAGGGCCcTGGACAC    AGG     2       CTLA4_site1:+500        Mismatches: 2
CTLA4_site9     CTLA4_site9_OB2_CTLA4_site2_500 GGAaTGAGGcCCATGGACAC    TGG     2       CTLA4_site2:+500        Mismatches: 2
CTLA4_site9     CTLA4_site9_OB2_CTLA4_site3_500 GGACTGgGGGCCtTGGACAC    AGG     2       CTLA4_site3:+500        Mismatches: 2
```
- The `Guide` column is the name of the on-target guide. This value is the same for all guides (because they all refer to the same on-target guide). Targets will be grouped by `Guide` value for plotting.
- The `Target` column specifies a customizable target name that will appear in reports and plots. Here, `MakeGuideFile` assigns each target a different name based on the Guide, the number of mismatches, and the chromosome and start position. In the super-small genome file, I named the chromosomes 'CTLA4_site0', 'CTLA4_site1', etc.
- The `Sequence` and `PAM` columns show the target sequence and PAM information.
- The `#MM` column shows the number of mismatches between the on-target and the off-target.
- The `Locus` column shows the genomic location of the target in the form Chromosome : Strand Position
- The `Mismatch_info` column is not required, but includes text describing the number of mismatches and bulges (if any). Note that Cas-offinder 3 is required for enumerating off-targets with bulges.

### Process
If you ran MakeGuideFile, you can now use the identified offtargets to run in Process mode using the command:
```
CRISPRessoSea Process --sample_file samples.demo.txt --target_file CRISPRessoSea_MakeGuideFileOutput/CRISPRessoSea.guide_info.txt --genome_file demo_genome.fa
```
- The `--sample_file` parameter specifies a text file specifying the names and fastq sequences of samples with columns `Name`, `fastq_r1`, and optionally `fastq_r2` and `group`.
- The `--target_file` parameter specifies a text file specifying the targets - in this case it is produced by the MakeGuideFile function.
- The `--genome_file` parameter specifies the path to the genome file.

If you didn't run that step, you can use the guide_info file in the demo dataset:
```
CRISPRessoSea Process --sample_file samples.demo.txt --target_file guides.demo.txt --genome_file demo_genome.fa
```

The `Process` function will produce an html report at CRISPRessoSea_output_on_samples.demo.txt/output_crispresso_sea.html, as well as an aggregated statistics file at CRISPRessoSea_output_on_samples.demo.txt/aggregated_stats_all.txt which can be modified and input for the Replot function below.
The aggregated_stats_all.txt file includes a row for each target. Columns specify information about each target (columns `target_id`, `target_name`, `target_chr`, `target_pos`, etc) as well as aggregated information from each sample (columns that end with `_highest_a_g_pct`, `_highest_c_t_pct`, `_highest_indel_pct` and `_tot_reads`). 

### Replot
If you'd like to change the order and name of guides, or add or change statistical tests, you can replot using the command:
```
CRISPRessoSea Replot --output_folder replot.output --reordered_stats_file replot_agg_stats.txt --reordered_sample_file replot_samples.txt --sig_method_parameters t_test,Control,Treated,0.05 
```
- The `--output_folder` parameter specifies where the replotted output will be produced.
- The `--reordered_stats_file` parameter is a modified aggregated stats file in the format produced by `Process`.
- The `--reordered_stats_file` parameter is a modified sample file where samples can be reordered if necessary.
- The `--sig_method_parameters` parameter specifies the significance test to be applied in the form of: 
 - none 
 - hard_cutoff,{cutoff}                     
 - mean_diff,{group1},{group2},{cutoff} 
 - t_test,{group1},{group2},{alpha} 
 - mann_whitney,{group1},{group2},{alpha}
 - neg_binomial,{group1},{group2},{alpha} 
 Here, we are using `t_test,Control,Treated,0.05` which specifies the t_test comparing Control vs Treated with a significance threshold of 0.05

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
usage: CRISPRessoSea.py Process [-h] [-o OUTPUT_FOLDER] -t TARGET_FILE -s SAMPLE_FILE [--gene_annotations GENE_ANNOTATIONS] -x GENOME_FILE [-r REGION_FILE] [-p N_PROCESSES]
                                [--crispresso_quantification_window_center CRISPRESSO_QUANTIFICATION_WINDOW_CENTER] [--crispresso_quantification_window_size CRISPRESSO_QUANTIFICATION_WINDOW_SIZE]
                                [--crispresso_base_editor_output] [--crispresso_default_min_aln_score CRISPRESSO_DEFAULT_MIN_ALN_SCORE] [--crispresso_plot_window_size CRISPRESSO_PLOT_WINDOW_SIZE]
                                [--allow_unplaced_chrs] [--plot_only_complete_targets] [--min_amplicon_coverage MIN_AMPLICON_COVERAGE] [--sort_based_on_mismatch]
                                [--allow_target_match_to_other_region_loc] [--top_percent_cutoff TOP_PERCENT_CUTOFF] [--min_amplicon_len MIN_AMPLICON_LEN] [--fail_on_pooled_fail]
                                [--plot_group_order PLOT_GROUP_ORDER] [--sig_method_parameters SIG_METHOD_PARAMETERS] [-v VERBOSITY] [--debug]

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
  -r REGION_FILE, --region_file REGION_FILE
                        Region file - a tab-separated .bed file with regions to analyze with columns for 'chr', 'start', and 'end'. If not provided, regions will be inferred by read alignment. (default:
                        None)
  -p N_PROCESSES, --n_processes N_PROCESSES
                        Number of processes to use. Set to "max" to use all available processors. (default: 8)
  --crispresso_quantification_window_center CRISPRESSO_QUANTIFICATION_WINDOW_CENTER
                        Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving
                        nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for
                        example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17 to only include mutations near the 5' end of the sgRNA. (default: -3)
  --crispresso_quantification_window_size CRISPRESSO_QUANTIFICATION_WINDOW_SIZE
                        Size (in bp) of the quantification window extending from the position specified by the '--cleavage_offset' or '--quantification_window_center' parameter in relation to the
                        provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the quantification window center are used in classifying reads as modified or unmodified. A value
                        of 0 disables this window and indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp. (default: 1)
  --crispresso_base_editor_output
                        Outputs plots and tables to aid in analysis of base editor studies. (default: False)
  --crispresso_default_min_aln_score CRISPRESSO_DEFAULT_MIN_ALN_SCORE
                        Default minimum homology score for a read to align to a reference amplicon. (default: 60)
  --crispresso_plot_window_size CRISPRESSO_PLOT_WINDOW_SIZE
                        Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the quantification_window_center for each guide are
                        plotted. (default: 20)
  --allow_unplaced_chrs
                        Allow regions on unplaced chromosomes (chrUn, random, etc). By default, regions on these chromosomes are excluded. If set, regions on these chromosomes will be included. (default:
                        False)
  --plot_only_complete_targets
                        Plot only targets with all values. If not set, all targets will be plotted. (default: False)
  --min_amplicon_coverage MIN_AMPLICON_COVERAGE
                        Minimum number of reads to cover a location for it to be plotted. Otherwise, it will be set as NA (default: 10)
  --sort_based_on_mismatch
                        Sort targets based on mismatch count. If true, the on-target will always be first (default: False)
  --allow_target_match_to_other_region_loc
                        If true, targets can match to regions even if the target chr:start is not in that region (e.g. if the target sequence is found in that region). If false/unset, targets can only
                        match to regions matching the target chr:start position. This flag should be set if the genome for guide design was not the same as the analysis genome. (default: False)
  --top_percent_cutoff TOP_PERCENT_CUTOFF
                        The top percent of aligned regions (by region read depth) to consider in finding non-overlapping regions during demultiplexing. This is a float between 0 and 1. For example, if
                        set to 0.2, the top 20% of regions (by read depth) will be considered. (default: 0.2)
  --min_amplicon_len MIN_AMPLICON_LEN
                        The minimum length of an amplicon to consider in finding non-overlapping regions during demultiplexing. Amplicons shorter than this will be ignored. (default: 50)
  --fail_on_pooled_fail
                        If true, fail if any pooled CRISPResso run fails. By default, processing will continue even if sub-CRISPResso commands fail. (default: False)
  --plot_group_order PLOT_GROUP_ORDER
                        Order of the groups to plot (if None, the groups will be sorted alphabetically) (default: None)
  --sig_method_parameters SIG_METHOD_PARAMETERS
                        Parameters for the significance method in the form of: none hard_cutoff,cutoff mean_diff,group1,group2,cutoff t_test,group1,group2,alpha mann_whitney,group1,group2,alpha
                        neg_binomial,group1,group2,alpha (default: None)
  -v VERBOSITY, --verbosity VERBOSITY
                        Verbosity level of output to the console (1-4) 4 is the most verbose (default: 3)
  --debug               Print debug information (default: False)
```

### Replot
```
usage: CRISPRessoSea.py Replot [-h] [-o OUTPUT_FOLDER] [-p FILE_PREFIX] -f REORDERED_STATS_FILE -s REORDERED_SAMPLE_FILE [--fig_width FIG_WIDTH] [--fig_height FIG_HEIGHT] [--seq_plot_ratio SEQ_PLOT_RATIO]
                               [--plot_group_order PLOT_GROUP_ORDER] [--sig_method_parameters SIG_METHOD_PARAMETERS] [--dot_plot_ylims DOT_PLOT_YLIMS] [--heatmap_max_value HEATMAP_MAX_VALUE]
                               [--heatmap_min_value HEATMAP_MIN_VALUE] [-v VERBOSITY] [--debug]

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
  --fig_width FIG_WIDTH
                        Width of the figure (default: 24)
  --fig_height FIG_HEIGHT
                        Height of the figure (default: 24)
  --seq_plot_ratio SEQ_PLOT_RATIO
                        Ratio of the width of the sequence plot to the data plot (>1 means the seq plot is larger than the data plot) (default: 1)
  --plot_group_order PLOT_GROUP_ORDER
                        Order of the groups to plot (if None, the groups will be sorted alphabetically) (default: None)
  --sig_method_parameters SIG_METHOD_PARAMETERS
                        Parameters for the significance method in the form of: none hard_cutoff,cutoff mean_diff,group1,group2,cutoff t_test,group1,group2,alpha mann_whitney,group1,group2,alpha
                        neg_binomial,group1,group2,alpha (default: None)
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
