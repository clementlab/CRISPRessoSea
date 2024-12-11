import argparse
import glob 
import numpy as np
import os
import pandas as pd
import seaborn as sns
import subprocess
import sys

from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPResso2Align

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def process_pools(sample_file, guide_file, genome_file, output_folder, n_processors=8, allow_unplaced_chrs=False, plot_only_complete_guides=False, 
                  min_amplicon_coverage=100, sort_based_on_mismatch=False, allow_guide_match_to_other_region_loc=False):
    """
    Discover amplicons and analyze sample editing at amplicons. This is the main entrypoint for this program.

    params:
    - sample_file: path to the sample file with headers: Name, fastq_r1, fastq_r2 (fastq_r2 is optional for single-end reads)
    - guide_file: path to the guide file with headers: Name, Sequence, PAM, Score, #MM, Gene, Locus
    - genome_file: path to the genome file (not including the trailing .fa)
    - output_folder: path to the output folder to produce data
    - n_processors: number of processors to use for multiprocessing
    - allow_unplaced_chrs: whether to allow regions to be identified on unplaced chromosomes (chrUn, random, etc). If true, these regions will be included.
    - plot_only_complete_guides: whether to plot only guides with data in all samples
    - min_amplicon_coverage: the minimum number of reads required to consider an amplicon
    - sort_based_on_mismatch: if true, guides are sorted based on mismatch count. If false, guides are presented in the order they are in the file
    - allow_guide_match_to_other_region_loc: if true, guides can be matched to regions based on sequence even if they are not in the same position. If false, guides can only match to regions containing that guide's chr:loc
    """
    if sample_file.endswith('.xlsx'):
        sample_df = pd.read_excel(sample_file)
    else:
        sample_df = pd.read_csv(sample_file, sep="\t")

    required_columns = ['Name','fastq_r1']

    for col in required_columns:
        if col not in sample_df.columns:
            raise Exception('Sample file must have column "' + col + '".\n' + 
                            'Found columns: ' + str(list(sample_df.columns)) + '\n' +
                            'Expecting columns ' + str(required_columns))
    
    #pull out first sample name which we will use to find amplicon locations
    first_sample = sample_df.iloc[0,:]
    first_sample_name = first_sample['Name']
    first_sample_r1 = first_sample['fastq_r1']
    first_sample_r2 = None
    if 'fastq_r2' in first_sample:
        first_sample_r2 = first_sample['fastq_r2']
    if genome_file.endswith('.fa'):
        genome_file = genome_file[:-3]
    
    for idx, row in sample_df.iterrows():
        sample_name = row['Name']
        sample_r1 = row['fastq_r1']
        sample_r2 = None
        if 'fastq_r2' in row and row['fastq_r2'] is not None:
            sample_r2 = row['fastq_r2']
        if not os.path.exists(sample_r1):
            raise Exception('Fastq R1 file for sample ' + sample_name + ' not found at ' + sample_r1)
        if sample_r2 is not None and not os.path.exists(sample_r2):
            raise Exception('Fastq R2 file for sample ' + sample_name + ' not found at ' + sample_r2)

    summary_output_folder = output_folder + 'summary/'
    if not os.path.exists(summary_output_folder):
        os.makedirs(summary_output_folder, exist_ok=True)

    crispresso_output_folder = output_folder + 'CRISPResso_output/'
    if not os.path.exists(crispresso_output_folder):
        os.makedirs(crispresso_output_folder, exist_ok=True)

    # first, run CRISPResso on the first sample to get the amplicon regions
    merged_regions, merged_regions_infos = run_initial_demux(first_sample_name, first_sample_r1, first_sample_r2, genome_file, output_folder=crispresso_output_folder, n_processors=n_processors, allow_unplaced_chrs=allow_unplaced_chrs)
    crispresso_region_file, guide_df, region_df = make_guide_region_assignments(merged_regions, merged_regions_infos, guide_file, output_folder, genome_file,
                                                                                sort_based_on_mismatch, allow_guide_match_to_other_region_loc=allow_guide_match_to_other_region_loc)

    aggregated_stats = guide_df.copy()
    aggregated_stats.set_index('guide_id', inplace=True)
    sample_df['CRISPRessoPooled_output_folder'] = 'NA'
    for sample_idx, sample_row in sample_df.iterrows():
        sample_name = sample_row['Name']
        sample_r1 = sample_row['fastq_r1']
        sample_r2 = None
        if 'fastq_r2' in sample_row:
            sample_r2 = sample_row['fastq_r2']
        this_pooled_run = run_crispresso_with_assigned_regions(sample_name, sample_r1, sample_r2, genome_file, crispresso_region_file, crispresso_output_folder, n_processors=n_processors)
        sample_df.loc[sample_idx, 'CRISPRessoPooled_output_folder'] = this_pooled_run
        if this_pooled_run is not None:
            guide_summary_file, guide_summary_df = analyze_run(summary_output_folder + sample_name, this_pooled_run, guide_df, region_df)
            sample_df.loc[sample_idx, 'guide_summary_file'] = guide_summary_file

            guide_summary_df.columns = ['guide_id','guide_label'] + [sample_name + '_' + x for x in guide_summary_df.columns[2:]]
            guide_summary_df.drop(columns=['guide_label'], inplace=True)

            aggregated_stats = pd.merge(aggregated_stats, guide_summary_df, how='left', on='guide_id')

            #for the columns in guide_summary_df, if tot_reads is less than min_amplicon_coverage, set all other columns to NA
            for col in guide_summary_df.columns:
                if col.startswith(sample_name) and col != sample_name + '_tot_reads':
                    aggregated_stats.loc[aggregated_stats[sample_name + '_tot_reads'] < min_amplicon_coverage, col] = np.nan

    aggregated_stats.sort_values(by='sort_index', inplace=True)
    aggregated_stats.to_csv(output_folder + 'aggregated_stats_all.txt', sep="\t", index=False)
    if plot_only_complete_guides:
        aggregated_stats_good = aggregated_stats.dropna()
    else:
        aggregated_stats_good = aggregated_stats
    guide_plot_df = create_guide_df_for_plotting(aggregated_stats_good)
    guide_plot_df.index = aggregated_stats_good['guide_chr'] + ':' + aggregated_stats_good['guide_pos'].astype(str) + ' ' + aggregated_stats_good['guide_name']

    print('Plotting for ' + str(len(aggregated_stats_good)) + '/' + str(len(aggregated_stats)) + ' guides after removing guides missing data in any samples')

    col_to_plot = 'highest_a_g_pct'
    col_title = 'Highest A/G %'
    plot_suffix = 'highest_a_g_pct.pdf'
    plot_guides_and_heatmap(guide_plot_df, aggregated_stats_good, col_to_plot, col_title, outfile_name=output_folder + plot_suffix)

    col_to_plot = 'highest_c_t_pct'
    col_title = 'Highest C/T %'
    plot_suffix = 'highest_c_t_pct.pdf'
    plot_guides_and_heatmap(guide_plot_df, aggregated_stats_good, col_to_plot, col_title, outfile_name=output_folder + plot_suffix)

    col_to_plot = 'highest_indel_pct'
    col_title = 'Highest Indel %'
    plot_suffix = 'highest_indel_pct.pdf'
    plot_guides_and_heatmap(guide_plot_df, aggregated_stats_good, col_to_plot, col_title, outfile_name=output_folder + plot_suffix)


def reverse_complement(seq):
    """
    Get the reverse complement of a sequence
    
    params:
    - seq: the sequence to reverse complement
    
    returns:
    - the reverse complement of the sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement[base] for base in reversed(seq))

def merge_locations(crispresso_pooled_genome_folder, genome_file, allow_unplaced_chrs=False, min_amplicon_len=50, debug=False):
    """
    Merge locations from a CRISPRessoPooled run to get a list of non-overlapping regions that could represent the original amplicons

    params:
    - crispresso_pooled_genome_folder: the folder containing the CRISPRessoPooled output
    - genome_file: path to the genome file
    - allow_unplaced_chrs: whether to allow regions on unplaced chromosomes (chrUn, random, etc). If true, these regions will be included.
    - min_amplicon_len: the minimum length of an amplicon to consider
    - debug: whether to print debug information

    returns:
    - list of merged regions
    - dictionary of merged regions -> region info (e.g. chr, start, end, read count)
    
    """
    genome_read_counts_file = os.path.join(crispresso_pooled_genome_folder, 'REPORT_READS_ALIGNED_TO_GENOME_ALL_DEPTHS.txt')
    if not os.path.isfile(genome_read_counts_file):
        raise Exception('Genome read counts file not found at ' + genome_read_counts_file)

    region_info = {}
    good_nonoverlapping_regions = []
    region_df = pd.read_csv(genome_read_counts_file, sep="\t")
    region_df_sub = region_df.sort_values(by='number of reads', ascending=False)
    top_percent_cutoff = 0.2
    top_percent_count = int(len(region_df_sub) * top_percent_cutoff)
    region_df_sub = region_df_sub.head(top_percent_count)
    print('Considering top ' + str(top_percent_cutoff * 100) + '% of regions (N=' + str(top_percent_count) + '/' + str(region_df.shape[0]) +') with at least ' + str(region_df_sub['number of reads'].min()) + \
        ' reads (mean reads is ' + str(region_df_sub['number of reads'].mean()) + ')')
    
    seen_region_seqs = {}
        
    for index, region_row in region_df_sub.iterrows():
        region_chr = region_row.chr_id
        region_start = region_row.start
        region_end = region_row.end
        region_count = region_row['number of reads']
        region_name = "_".join([region_chr, str(region_start), str(region_end)])
        region_seq_output = subprocess.check_output(
            '%s faidx %s %s:%d-%d'%('samtools',genome_file + ".fa",region_chr,region_start,region_end-1),shell=True).decode(sys.stdout.encoding)
        region_seq = "".join(region_seq_output.split("\n")[1:])

        region_info[region_name] = {'chr': region_chr, 'start': region_start, 'end': region_end, 'region_count': region_count, 'seq': region_seq}
        
        if not allow_unplaced_chrs:
            if 'chrUn' in region_chr or 'random' in region_chr:
                continue
        region_len = region_end - region_start
        if region_len < min_amplicon_len:
            continue

        if region_seq.upper() in seen_region_seqs or reverse_complement(region_seq.upper()) in seen_region_seqs:
            continue

        overlapped_any_region = False
        for other_region in good_nonoverlapping_regions:
            other_region_chr = region_info[other_region]['chr']
            other_region_start = region_info[other_region]['start']
            other_region_end = region_info[other_region]['end']
            overlaps_this_region = False
            if region_chr == other_region_chr:
                if other_region_start <= region_start <= other_region_end:
                    overlaps_this_region = True
                elif other_region_start <= region_end <= other_region_end:
                    overlaps_this_region = True
                elif other_region_start >= region_start and other_region_end <= region_end:
                    overlaps_this_region = True
                elif other_region_start <= region_start and other_region_end >= region_end:
                    overlaps_this_region = True
            if overlaps_this_region:
                overlapped_any_region = True
                other_region_count = region_info[other_region]['region_count']
                if region_count > other_region_count:
                    good_nonoverlapping_regions.remove(other_region)
                    good_nonoverlapping_regions.append(region_name)
        if not overlapped_any_region:
            good_nonoverlapping_regions.append(region_name)
            seen_region_seqs[region_seq.upper()] = True

    total_region_count = region_df.shape[0]
    kept_region_count = len(good_nonoverlapping_regions)
    print('Kept ' + str(kept_region_count) + '/' + str(total_region_count) + ' regions')
    if debug:
        for region in good_nonoverlapping_regions:
            print(region)
    good_nonoverlapping_region_infos = {}
    for region in good_nonoverlapping_regions:
        this_region_info = region_info[region]
        good_nonoverlapping_region_infos[region] = this_region_info
    return good_nonoverlapping_regions, good_nonoverlapping_region_infos

def run_initial_demux(experiment_name, fastq_r1, fastq_r2, genome_file, output_folder='', n_processors=8, allow_unplaced_chrs=False, suppress_output=True):
    """
    Run CRISPResso on input reads to determine which regions are frequently aligned to

    params:
    - experiment_name: a name for the experiment
    - fastq_r1: path to the R1 fastq file
    - fastq_r2: path to the R2 fastq file
    - genome_file: path to the genome file
    - output_folder: path to the output folder where results should be written
    - n_processors: number of processors to use
    - allow_unplaced_chrs: whether to allow regions on bad chromosomes (chrUn, random, etc). If true, these regions will be included.

    returns:
    - list of merged regions
    - dictionary of merged regions -> region info (e.g. chr, start, end, read count)
    """

    r2_string = ''
    if fastq_r2 is not None:
        r2_string = ' -r2 ' + fastq_r2
    output_string = ''
    if output_folder != '':
        output_string = ' -o ' + output_folder
    suppress_output_string = ''
    if suppress_output:
        suppress_output_string = ' --verbosity 0'

    CRISPResso_output_folder = output_folder + "CRISPRessoPooled_on_"+experiment_name+"_demux"
    command = 'CRISPRessoPooled -x ' + genome_file + ' -r1 ' + fastq_r1 + r2_string + output_string + ' -n ' + experiment_name + '_demux -p ' + str(n_processors) + ' --no_rerun --keep_intermediate --suppress_plots' + suppress_output_string

    crispresso_run_is_complete = False
    if os.path.exists(CRISPResso_output_folder):
        try:
            crispresso_pooled_info = CRISPRessoShared.load_crispresso_info(crispresso_info_file_path=CRISPResso_output_folder + '/CRISPResso2Pooled_info.json')
            if 'demultiplexing_genome_only_regions' in crispresso_pooled_info['running_info']['finished_steps']:
                crispresso_run_is_complete = True
                print('CRISPResso output folder already exists for initial demultiplexing, skipping CRISPResso run')
        except:
            pass
    if not crispresso_run_is_complete:
        print('Aligning reads to the genome to find amplicon locations')
        print('Running command ' + str(command))
        subprocess.run(command, shell=True, check=True)

    merged_regions, merged_regions_infos = merge_locations(CRISPResso_output_folder, genome_file, allow_unplaced_chrs=allow_unplaced_chrs)

    return merged_regions, merged_regions_infos

def parse_guide_info(guide_file, sort_based_on_mismatch=False):
    """
    Parse a guide file in tab-separated or xls format. Required columns are:
    Name: the name of the on-target guide
    Sequence: the sequence of the target site (either on-target or off-target)
    PAM: the PAM at the target site
    Score: the score of the off-target
    #MM: the number of mismatches between the on- and off-target. If this is 0 or NA, this guide is considered the on-target
    Locus: the locus of the target site in the genome, in the format chr:+pos for positive strand or chr:-pos for negative strand 
    
    Optional column:
    anno: an annotation name for the guide

    params:
    guide_file: path to the guide file
    sort_based_on_mismatch: if true, guides are sorted based on mismatch count. If false, guides are presented in the order they are in the file

    The On-target for each guide is the guide with the fewest mismatches (may not necessarily be 0MM for retargeting guides). If there are multiple guides with the same number of mismatches, the first one is used.

    Returns:
    - a pandas dataframe with the following columns for processing:
        guide_id: a unique identifier for the guide
        guide_name: the name of the guide
        sort_index: the index of the guide in the original file
        guide_chr: the chromosome of the guide
        guide_pos: the position of the guide
        guide_seq_no_gaps_with_pam: the guide sequence without gaps and with PAM
        guide_seq_no_gaps: the guide sequence without gaps
        guide_seq_with_gaps: the guide sequence with gaps
        guide_pam: the guide pam sequence
        ontarget_name: the name of the on-target guide
        ontarget_sequence: the sequence of the on-target guide

    - a dictionary of guide_id -> guide_name

    """
    guide_name_lookup = {}

    if not os.path.exists(guide_file):
        raise Exception('Guide file not found at ' + guide_file)

    if guide_file.endswith('.xlsx'):
        guide_df = pd.read_excel(guide_file)
    else:
        guide_df = pd.read_csv(guide_file, sep="\t")

    required_columns = ['Name','Sequence','PAM','Score','#MM','Gene','Locus']

    for col in required_columns:
        if col not in guide_df.columns:
            raise Exception('Guide file must have column "' + col + '".\n' + 
                            'Found columns: ' + str(list(guide_df.columns)) + '\n' +
                            'Expecting columns ' + str(required_columns))
    
    guide_df['#MM'] = guide_df['#MM'].fillna(0).astype(int)
    guide_df['sort_index'] = guide_df.index

    on_targets = guide_df.sort_values(by=['#MM', 'sort_index']).drop_duplicates(subset='Name', keep='first')
    if sort_based_on_mismatch:
        guide_df['sort_index'] = list(range(1,guide_df.shape[0]))

    for idx, row in on_targets.iterrows():
        if row['#MM'] != 0:
            print('Warning: On-target guide ' + row['Name'] + ' has ' + str(row['#MM']) + ' mismatches')    

    ontarget_seqs = {}
    for idx, row in on_targets.iterrows():
        ontarget_seqs[row['Name']] = row['Sequence']

    guide_df['ontarget_name'] = 'NA'
    guide_df['ontarget_sequence'] = 'NA'
    guide_df['guide_id'] = 'NA' #unique id
    guide_df['guide_name'] = 'NA' #may be supplied by user
    guide_df['guide_chr'] = 'NA'
    guide_df['guide_pos'] = 'NA'
    guide_df['guide_seq_with_gaps'] = guide_df['Sequence']
    guide_df['guide_pam'] = guide_df['PAM']
    guide_df['guide_seq_no_gaps'] = 'NA'
    guide_df['guide_seq_no_gaps_with_pam'] = 'NA'
    for idx, row in guide_df.iterrows():
        this_ontarget_name = row['Name']
        if this_ontarget_name not in ontarget_seqs:
            raise Exception('On-target sequence not found for guide ' + this_ontarget_name + '. Found ontargets: ' + str(list(ontarget_seqs.keys())))
        if str(row['#MM']) == '0':
            this_guide_id = this_ontarget_name + "_ON"
        else:
            this_guide_id = this_ontarget_name + "_OB" + str(int(row['#MM']))

        guide_df.loc[idx,'guide_id'] = str(idx) + ' ' + this_guide_id

        this_guide_name = this_guide_id
        if 'anno' in row and row['anno'] is not None and str(row['anno']) != 'nan': 
            this_guide_name = row['anno']
        guide_df.loc[idx,'guide_name'] = this_guide_name

        if this_ontarget_name not in ontarget_seqs:
            raise Exception('Cannot find ontarget name for ' + str(this_ontarget_name) + ' in ' + str(ontarget_seqs.keys()))
        this_ontarget_seq = ontarget_seqs[this_ontarget_name]
        guide_df.loc[idx,'ontarget_sequence'] = this_ontarget_seq
        guide_df.loc[idx,'ontarget_name'] = this_ontarget_name
        
        this_chr_loc_els = row['Locus'].split(':')
        this_chr = this_chr_loc_els[0]
        this_pos = int(this_chr_loc_els[1][1:])
        guide_df.loc[idx,'guide_chr'] = this_chr
        guide_df.loc[idx,'guide_pos'] = this_pos
        guide_df.loc[idx,'guide_seq_no_gaps_with_pam'] = row['Sequence'].replace('-','') + row['PAM']
        guide_df.loc[idx,'guide_seq_no_gaps'] = row['Sequence'].replace('-','')

        guide_name_lookup[this_guide_id] = this_guide_name
    
    return guide_df[['guide_id','guide_name','sort_index','guide_chr','guide_pos','guide_seq_no_gaps_with_pam','guide_seq_no_gaps','guide_seq_with_gaps', 'guide_pam', 'ontarget_name', 'ontarget_sequence']], guide_name_lookup



def make_guide_region_assignments(merged_regions, merged_regions_infos, guide_file, output_folder, genome_file, sort_based_on_mismatch=False, allow_guide_match_to_other_region_loc=False):
    """
    Assign guides to regions based on their sequences and positions

    params:
    - merged_regions: list of merged regions
    - merged_regions_infos: dictionary of merged regions -> region info (e.g. chr, start, end, read count)
    - guide_file: path to the user-provided guide file
    - output_folder: path to the output file root
    - genome_file: path to the genome file to look up region seqeunces
    - sort_based_on_mismatch: if true, guides are sorted based on mismatch count. If false, guides are presented in the order they are in the file
    - allow_guide_match_to_other_region_loc: if true, guides can be matched to regions based on sequence even if they are not in the same position. If false, guides can only match to regions containing that guide's chr:loc

    returns:
    - path to the CRISPRessoPooledRegions file 
    - guide_df: the guide dataframe
    - region_df: the region dataframe
    """
    #we parse the guide_df here because we'll add more columns to it
    guide_df, guide_name_lookup = parse_guide_info(guide_file, sort_based_on_mismatch)

    for region in merged_regions:
        region_chr = merged_regions_infos[region]['chr']
        region_start = merged_regions_infos[region]['start']
        region_end = merged_regions_infos[region]['end']
        region_seq = merged_regions_infos[region]['seq']

    # first, for each guide, match its sequence to genomic region(s)
    # the matching will be performed so that each guide is assigned to
    # 1) If there is one region that corresponds to the guide position it is assigned to that region.
    # 2) Else if there are multiple regions that correspond to the guide position, it is assigned to the region with the highest read count
    # 3) Else (the chr:pos is not given or the guide position is in no regions), if the guide sequence is found in one region the guide is assigned to that region
    # 4) Else (there are multiple regions that contain the guide sequence), it is assigned to the region with the highest read count without a previously-assigned guide.
    guide_matches = {} # guide_seq_id -> final assigned region
    region_matches = defaultdict(list) # region -> list of guide_seq_ids that were assigned to it

    all_guide_matches = defaultdict(list) # guide_seq_id -> list of regions (all regions with position or sequence match)
    all_region_matches = defaultdict(list) # region -> list of guide_seq_ids

    for guide_idx, guide_row in guide_df.iterrows():
        guide_seq_with_pam = guide_row['guide_seq_no_gaps_with_pam']
        guide_id = guide_row['guide_id']
        guide_name = guide_row['guide_name']
        guide_chr = guide_row['guide_chr']
        guide_pos = guide_row['guide_pos']

        this_guide_seq_matches = []
        this_guide_region_matches = []
        for region in merged_regions:
            region_info = merged_regions_infos[region]
            region_chr = region_info['chr']
            region_start = region_info['start']
            region_end = region_info['end']
            region_seq = region_info['seq']
            
            # first, check if the region contains the guide sequence
            is_seq_match = False
            if guide_seq_with_pam.upper() in region_seq.upper():
                is_seq_match = True
            elif reverse_complement(guide_seq_with_pam).upper() in region_seq.upper():
                is_seq_match = True

            if not allow_guide_match_to_other_region_loc: #don't allow matches based on guide match to region sequence
                is_seq_match = False

            if is_seq_match:
                this_guide_seq_matches.append(region)
            
            # next, check if the guide position is within the region
            is_pos_match = False
            if guide_chr == region_chr and guide_pos >= region_start and guide_pos <= region_end:
                this_guide_region_matches.append(region) 
                is_pos_match = True
            elif guide_chr.replace('chr','') == region_chr.replace('chr','') and guide_pos >= region_start and guide_pos <= region_end:
                this_guide_region_matches.append(region) 
                is_pos_match = True

            if is_seq_match or is_pos_match:
                all_guide_matches[guide_id].append(region)
                all_region_matches[region].append(guide_id)
            
            
        # make the final assignment for this guide
        if len(this_guide_region_matches) == 1: # if there is one region match, take it
            guide_matches[guide_id] = this_guide_region_matches[0]
            region_matches[this_guide_region_matches[0]].append(guide_id)
        elif len(this_guide_region_matches) > 1: # if there are multiple region matches, take the one with highest number of reads
            best_region = None
            best_region_count = 0
            for region in this_guide_region_matches:
                region_count = merged_regions_infos[region]['region_count']
                if region_count > best_region_count:
                    best_region = region
                    best_region_count = region_count
            guide_matches[guide_id] = best_region
            region_matches[best_region].append(guide_id)
        elif len(this_guide_seq_matches) == 1: # if there is one sequence match, take it
            guide_matches[guide_id] = this_guide_seq_matches[0]
            region_matches[this_guide_seq_matches[0]].append(guide_id)
        elif len(this_guide_seq_matches) > 1: # if there are multiple sequence matches, take the one with highest number of reads without a guide match
            best_region = this_guide_seq_matches[0]
            best_region_count = 0
            for region in this_guide_seq_matches:
                region_count = merged_regions_infos[region]['region_count']
                if region_count > best_region_count and len(region_matches[region]) == 0:
                    best_region = region
                    best_region_count = region_count
            guide_matches[guide_id] = best_region
            region_matches[best_region].append(guide_id)

    # next, finalize matches and write them to file
    out_file = output_folder + 'guide_matches.txt'
    guide_match_count = 0 # how many guides had a region match 
    guide_nomatch_count = 0  # how many guides did not have a region match

    guide_df['matched_region_count'] = 0
    guide_df['matched_region'] = 'NA'
    with open(out_file, 'w') as fout:
        fout.write('guide_id\tguide_name\tguide_seq\tguide_seq_with_pam\tguide_chr\tguide_pos\tontarget_name\tontarget_seq\tmatched_region_count\tmatched_regions\tregion_chr\tregion_start\tregion_end\tregion_read_count\tregion_seq\n')
        for guide_idx, guide_row in guide_df.iterrows():
            guide_seq_id = guide_row['guide_id']
            guide_seq = guide_row['guide_seq_no_gaps']
            guide_seq_with_pam = guide_row['guide_seq_no_gaps_with_pam']
            guide_seq_name = guide_row['guide_name']
            guide_chr = guide_row['guide_chr']
            guide_pos = guide_row['guide_pos']
            guide_ontarget_name = guide_row['ontarget_name']
            guide_ontarget_seq = guide_row['ontarget_sequence']
            if guide_seq_id in guide_matches:
                guide_match_count += 1
                matched_region_count = len(all_guide_matches[guide_seq_id])
                matched_regions = ', '.join(all_guide_matches[guide_seq_id])
                region_name = guide_matches[guide_seq_id]
                region_info = merged_regions_infos[region_name]
                region_chr = region_info['chr']
                region_start = region_info['start']
                region_end = region_info['end']
                region_seq = region_info['seq']
                fout.write('\t'.join([str(x) for x in [guide_seq_id, guide_seq_name, guide_seq, guide_seq_with_pam, guide_chr, guide_pos, guide_ontarget_name, guide_ontarget_seq, matched_region_count, matched_regions, region_chr, region_start, region_end, region_info['region_count'], region_seq]]) + '\n')
                guide_df.loc[guide_idx, 'matched_region_count'] = matched_region_count
                guide_df.loc[guide_idx, 'matched_region'] = region_name
            # if guide wasn't selected as the final match for a region, it will be in all_guide_matches
            elif guide_seq_id in all_guide_matches:
                guide_nomatch_count += 1
                matched_region_count = len(all_guide_matches[guide_seq_id])
                matched_regions = ', '.join(all_guide_matches[guide_seq_id])
                fout.write('\t'.join([str(x) for x in [guide_seq_id, guide_seq_name, guide_seq, guide_seq_with_pam, guide_chr, guide_pos, guide_ontarget_name, guide_ontarget_seq, 0, matched_regions, 'NA', 'NA', 'NA', 'NA', 'NA']]) + '\n')
            else:
                guide_nomatch_count += 1
                fout.write('\t'.join([str(x) for x in [guide_seq_id, guide_seq_name, guide_seq, guide_seq_with_pam, guide_chr, guide_pos, guide_ontarget_name, guide_ontarget_seq, 0, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']]) + '\n')

    print('Matched ' + str(guide_match_count) + '/' + str(len(guide_df)) + ' guides to regions from read alignments. Wrote matches to ' + out_file)

    # next write information about all regions plus the file for crispresso input
    crispresso_output_file = output_folder + 'CRISPRessoPooledRegions.txt'
    all_region_output_file = output_folder + 'all_regions.txt'
    region_match_count = 0
    region_nomatch_count = 0
    printed_crispresso_count = 0
    with open(all_region_output_file,'w') as rout, open(crispresso_output_file,'w') as cout:
        rout.write('region_id\tchr\tstart\tend\tread_count\tseq\tguide_match_count\tguide_matches\tguide_id\tguide_name\tguide_seq\tguide_chr\tguide_pos\n')
        for region in merged_regions:
            region_info = merged_regions_infos[region]
            region_chr = region_info['chr']
            region_start = region_info['start']
            region_end = region_info['end']
            region_seq = region_info['seq']
            region_name = "_".join([region_chr, str(region_start), str(region_end)])
            guide_id = 'NA'
            guide_name = 'NA'
            guide_seq = 'NA'
            guide_chr = 'NA'
            guide_pos = 'NA'
            this_guide_match_count = len(all_region_matches[region])
            if this_guide_match_count > 0:
                this_guide_matches = ', '.join(all_region_matches[region])
            else:
                this_guide_matches = 'NA'
            if region in all_region_matches:
                if region in region_matches:
                    region_match_count += 1
                    guide_id = region_matches[region][0] # the assigned guides that matched to this region
                    guide_seq = guide_df.loc[guide_df['guide_id'] == guide_id, 'guide_seq_no_gaps'].values[0] # no pam (for CRISPResso)
                    guide_name = guide_df.loc[guide_df['guide_id'] == guide_id, 'guide_name'].values[0]
                    guide_chr = guide_df.loc[guide_df['guide_id'] == guide_id, 'guide_chr'].values[0]
                    guide_pos = guide_df.loc[guide_df['guide_id'] == guide_id, 'guide_pos'].values[0]
                    region_name = region_name + "_" + guide_name
                    cout.write("\t".join([region_name,region_seq,guide_seq]) + '\n')
                    printed_crispresso_count += 1
                else:
                    region_nomatch_count += 1
            else:
                region_nomatch_count += 1
            rout.write('\t'.join(str(x) for x in [region_name, region_chr, region_start, region_end, region_info['region_count'], region_seq, this_guide_match_count, this_guide_matches, guide_id, guide_name, guide_seq, guide_chr, guide_pos]) + '\n')
        
    print('Matched ' + str(region_match_count) + "/" + str(region_match_count + region_nomatch_count) + ' frequently-aligned locations to guides. Wrote region info to ' + all_region_output_file)

    if printed_crispresso_count == 0:
        raise Exception('No regions matched to guides. Check the guide file and the regions file for consistency')
    region_df = pd.read_csv(all_region_output_file, sep="\t")
    return crispresso_output_file, guide_df, region_df


def run_crispresso_with_assigned_regions(experiment_name, fastq_r1, fastq_r2, genome_file, crispresso_region_file, output_folder, n_processors=8, suppress_output=True, skip_failed_subruns=True):
    """
    Run CRISPResso on a specific sample with the assigned regions

    params:
    - experiment_name: a name for the experiment
    - fastq_r1: path to the R1 fastq file
    - fastq_r2: path to the R2 fastq file
    - genome_file: path to the genome file
    - crispresso_region_file: path to the file with region assignments
    - n_processors: number of processors to use
    - suppress_output: whether to suppress output from CRISPRessoPooled
    - skip_failed_subruns: whether to skip failed CRISPRessoPooled subruns. If true, CRISPRessoPooled will complete even if individual subruns fail (e.g. due to too few aligned reads).

    returns:
    - the output folder for this CRISPResso run
    """
    # Run CRISPResso again on the subset, this time including guide sequences
    r2_string = ''
    if fastq_r2 is not None:
        r2_string = ' -r2 ' + fastq_r2
    output_string = ''
    if output_folder != '':
        output_string = ' -o ' + output_folder

    suppress_output_string = ''
    if suppress_output:
        suppress_output_string = ' --verbosity 0'

    skip_failed_string = ""
    if skip_failed_subruns:
        skip_failed_string = " --skip_failed"

    CRISPResso_output_folder = output_folder + "CRISPRessoPooled_on_"+experiment_name

    command = 'CRISPRessoPooled -f ' + crispresso_region_file + ' -x ' + genome_file + ' -n ' + experiment_name + ' ' + ' -r1 ' + fastq_r1 + r2_string + output_string + ' --min_reads_to_use_region 1 --default_min_aln_score 20 ' + \
                ' --base_editor_output --quantification_window_center -10 --quantification_window_size 16' + \
                ' -p ' + str(n_processors) + ' --no_rerun  --exclude_bp_from_left 0 --exclude_bp_from_right 0 --plot_window_size 10' + suppress_output_string + skip_failed_string

    crispresso_run_is_complete = False
    if os.path.exists(CRISPResso_output_folder):
        try:
            crispresso_pooled_info = CRISPRessoShared.load_crispresso_info(crispresso_info_file_path=CRISPResso_output_folder + '/CRISPResso2Pooled_info.json')
            if 'end_time' in crispresso_pooled_info['running_info']:
                crispresso_run_is_complete = True
                print('CRISPResso output already exists for sample ' + experiment_name + ', skipping CRISPResso run')
        except:
            raise Exception('failed!')
            pass

    if not crispresso_run_is_complete:
        print('Aligning reads to the genome to find amplicon locations')
        print('Running command ' + str(command))
        try:
            subprocess.run(command, shell=True, check=True)
        except:
            print('Failed to run CRISPRessoPooled for sample ' + experiment_name)
            return None

    return CRISPResso_output_folder


def analyze_run(output_name, crispresso_folder, guide_df, region_df):
    """
    Analyzes a CRISPResso batch run to get the base editing and indel rates for each site

    params:
    - output_name: the name of the output file
    - crispresso_folder: the folder containing the CRISPResso output
    - guide_df: the guide dataframe
    - region_df: the region dataframe

    returns: 
    - summary_output_file: a file with the summary of the run
    - guide_summary_df: a pandas data frame with results for each guide
    """

    guide_summary_file = output_name + ".complete_guide_summary.txt"
    if crispresso_folder is None:
        guide_data = []
        with open(guide_summary_file, 'w') as fout:
            fout.write('guide_id\tguide_label\tpooled_result_name\thighest_a_g_pct\thighest_c_t_pct\thighest_indel_pct\ttot_reads\n')
        for idx, row in guide_df.iterrows():
            guide_id = row['guide_id']
            guide_name = row['guide_name']
            fout.write(guide_id + '\t' + guide_name + '\tNA\tNA\tNA\tNA\n')
            guide_data.append([guide_id, guide_name, np.nan, np.nan, np.nan, np.nan])
        guide_data_df = pd.DataFrame(guide_data, columns=['guide_id', 'guide_label', 'pooled_result_name', 'highest_a_g_pct', 'highest_c_t_pct', 'highest_indel_pct', 'tot_reads'])
        return guide_summary_file, guide_data_df

    pooled_info = CRISPRessoShared.load_crispresso_info(crispresso_info_file_path=crispresso_folder + '/CRISPResso2Pooled_info.json')
    pooled_results = pooled_info['results']['good_region_folders']

    output_dir = output_name + "_plots"
    os.makedirs(output_dir, exist_ok=True)

    guide_results = {}
    output_summary_file = output_name + ".complete_summary.txt"
    with open(output_summary_file, 'w') as fout:
        fout.write("\t".join(['folder_name', 'highest_a_g_pct', 'highest_c_t_pct', 'highest_indel_pct', 'tot_reads', 'guide_id','guide_name', 'guide_chr','guide_pos', 'guide_used', 'guide_seq_with_gaps', 'guide_pam']) + '\n')
        for pooled_result_name in pooled_results.keys():
            folder_name = pooled_results[pooled_result_name]
            run_info = CRISPRessoShared.load_crispresso_info(crispresso_info_file_path=crispresso_folder + "/" + folder_name + '/CRISPResso2_info.json')
            guide_used = run_info['results']['refs']['Reference']['sgRNA_sequences'][0]

            region_info = region_df.loc[region_df['region_id'] == pooled_result_name]
            guide_id = region_info['guide_id'].values[0]
            guide_info = guide_df.loc[guide_df['guide_id'] == guide_id]
            guide_name = guide_info['guide_name'].values[0]

            include_idxs = run_info['results']['refs']['Reference']['include_idxs']
            nuc_freq_table_file = crispresso_folder + "/" + folder_name + '/Nucleotide_frequency_table.txt'
            mod_freq_table_file = crispresso_folder + "/" + folder_name + '/Modification_count_vectors.txt'
            nuc_freq_table = pd.read_csv(nuc_freq_table_file, sep='\t', index_col=0)
            mod_freq_table = pd.read_csv(mod_freq_table_file, sep='\t', index_col=0)
            highest_a_g_pct = 0
            highest_c_t_pct = 0
            highest_indel_pct = 0
            this_tot = 0
            for this_idx in include_idxs:
#        for i in range(8, 24): # starts and ends 3bp from the end of the guide (skips subs in the 2bp on the very ends)
#                this_idx = include_idxs[i]
                nuc_vals = nuc_freq_table.iloc[:,this_idx]
                mod_vals = mod_freq_table.iloc[:,this_idx]

                this_nuc = nuc_vals.name
                a_count = nuc_vals['A']
                c_count = nuc_vals['C']
                g_count = nuc_vals['G']
                t_count = nuc_vals['T']
                del_count = mod_vals['Deletions']
                ins_count = mod_vals['Insertions']
                tot_count = mod_vals['Total']
                this_tot = tot_count

                if this_nuc[0] == 'A':
                    a_g_pct = g_count / tot_count
                    if a_g_pct > highest_a_g_pct:
                        highest_a_g_pct = a_g_pct
                elif this_nuc[0] == 'C':
                    c_t_pct = t_count / tot_count
                    if c_t_pct > highest_c_t_pct:
                        highest_c_t_pct = c_t_pct
                elif this_nuc[0] == 'G':
                    c_t_pct = a_count / tot_count
                    if c_t_pct > highest_c_t_pct:
                        highest_c_t_pct = c_t_pct
                elif this_nuc[0] == 'T':
                    a_g_pct = c_count / tot_count
                    if a_g_pct > highest_a_g_pct:
                        highest_a_g_pct = a_g_pct
                if (del_count+ins_count)/tot_count > highest_indel_pct:
                    highest_indel_pct = (del_count+ins_count)/tot_count
            highest_a_g_pct *= 100
            highest_c_t_pct *= 100
            highest_indel_pct *= 100
            fout.write("\t".join([str(x) for x in [folder_name, highest_a_g_pct, highest_c_t_pct, highest_indel_pct, this_tot, guide_id, 
                                                   guide_info['guide_name'].values[0], guide_info['guide_chr'].values[0], guide_info['guide_pos'].values[0], 
                                                   guide_used, guide_info['guide_seq_with_gaps'].values[0], guide_info['guide_pam'].values[0]]]) + '\n')
            guide_results[guide_id] = (pooled_result_name, highest_a_g_pct, highest_c_t_pct, highest_indel_pct, this_tot)

            plot2b_file = crispresso_folder + "/" + folder_name + '/' + run_info['results']['refs']['Reference']['plot_2b_roots'][0]+".pdf"
            if not os.path.isfile(plot2b_file):
                print('WARNING: missing plot2b file for', pooled_result_name)
                continue
            destination_name = guide_name + "." + pooled_result_name + ".plot2b.pdf"
            os.popen('cp ' + plot2b_file + ' ' + output_dir + "/" + destination_name) 

    guide_data = []
    with open(guide_summary_file, 'w') as fout:
        fout.write('guide_id\tguide_label\tpooled_result_name\thighest_a_g_pct\thighest_c_t_pct\thighest_indel_pct\ttot_reads\n')
        for idx, row in guide_df.iterrows():
            guide_id = row['guide_id']
            guide_name = row['guide_name']
            if guide_id not in guide_results:
                fout.write(guide_id + '\t' + guide_name + '\tNA\tNA\tNA\tNA\n')
                guide_data.append([guide_id, guide_name, np.nan, np.nan, np.nan, np.nan])
            else:
                pooled_result_name, highest_a_g_pct, highest_c_t_pct, highest_indel_pct, tot_reads = guide_results[guide_id]
                fout.write(guide_id + '\t' + guide_name + '\t' + pooled_result_name + '\t' + str(highest_a_g_pct) + '\t' + str(highest_c_t_pct) + '\t' + str(highest_indel_pct) + '\t' + str(tot_reads) + '\n')
                guide_data.append([guide_id, guide_name, pooled_result_name, highest_a_g_pct, highest_c_t_pct, highest_indel_pct, tot_reads])
    plot_heatmap(output_name)
    guide_data_df = pd.DataFrame(guide_data, columns=['guide_id', 'guide_label', 'pooled_result_name', 'highest_a_g_pct', 'highest_c_t_pct', 'highest_indel_pct', 'tot_reads'])
    return guide_summary_file, guide_data_df

def view_complete_guide_summary(output_name):
    """
    View the complete_guide_summary.txt file for debugging

    params:
    - output_name: the name of the output folder

    """
    d = pd.read_csv(output_name + ".complete_guide_summary.txt", sep='\t')
    d.dropna(inplace=True)
    d.index = d['guide_id']
    print(d)

def plot_heatmap(output_name):
    """
    Plot a heatmap of the highest_a_g_pct, highest_c_t_pct, and highest_indel_pct for each guide
    Data is read from output_name + ".complete_guide_summary.txt"
    Plot is written to output_name + ".heatmap.pdf"

    params:
    - output_name: the name of the output folder
    """
    d = pd.read_csv(output_name + ".complete_guide_summary.txt", sep='\t')
    d.dropna(inplace=True)
    d.set_index('guide_id', inplace=True)
    fig = plt.figure(figsize=(12,12), dpi= 100, facecolor='w', edgecolor='k')
    sns.heatmap(d[['highest_a_g_pct', 'highest_c_t_pct', 'highest_indel_pct']].astype(float), annot=True, fmt=".4f",cmap='Reds')
    fig.savefig(output_name + ".heatmap.pdf") 
    plt.close(fig)

def plot_heatmap_sub(output_name, guide_name):
    """
    Plot a heatmap of the highest_a_g_pct, highest_c_t_pct, and highest_indel_pct for a specific guide
    Data is read from output_name + ".complete_guide_summary.txt"
    Plot is written to output_name + ".heatmap.pdf"

    params:
    - output_name: the name of the output folder
    - guide_name: the name of the guide to plot
    """
    d = pd.read_csv(output_name + ".complete_guide_summary.txt", sep='\t')
    d.dropna(inplace=True)
    d.set_index('guide_id', inplace=True)
    dsub = d[d['ontarget_name'] == guide_name]
    fig = plt.figure(figsize=(6,6), dpi= 100, facecolor='w', edgecolor='k')
    sns.heatmap(dsub[['highest_a_g_pct', 'highest_c_t_pct', 'highest_indel_pct']].astype(float), annot=True, fmt=".4f",cmap='Reds')
    plt.close(fig)

def create_guide_df_for_plotting(guide_df):
    """
    Create a dataframe for plotting guides in a heatmap
    
    params:
    - guide_df: the guide dataframe
        includes columns: guide_id, guide_name, guide_seq_no_gaps, guide_pam, ontarget_name, ontarget_sequence
    
    returns:
    - a dataframe with guide sequences for plotting
        adds column: padded_seqs_for_plot (list of characters for each box in the plot heatmap)
    """
    aln_matrix = CRISPResso2Align.make_matrix()

    names_for_plot = []
    seqs_for_plot = []
    aln_ontarget_seqs = []
    aln_offtarget_seqs = []

    last_guide_name = None
    last_guide_seq = None
    for idx, row in guide_df.iterrows():
        if row['ontarget_name'] == last_guide_name:
            ref_incentive = np.zeros(len(last_guide_seq) + 1, dtype=int)
            off_target_aln, on_target_aln, aln_score = CRISPResso2Align.global_align(row['guide_seq_no_gaps'], last_guide_seq, matrix=aln_matrix,
                                                                        gap_incentive=ref_incentive,
                                                                        gap_open=-15,
                                                                        gap_extend=-10, )
            
            this_chars_for_plot = []
            first_base_is_del_in_ontarget = False # toggle if starts with run of deletions in ontarget
            # bases are collapsed left meaning that an insertion will be grouped with the base to its left (except for the case when a deletion starts, and it is grouped to the right)
            # example:
            #  on_target_aln  = "---A-A-CCTG--"
            #  off_target_aln = "ACTAGATCCTGCC"
            # produces:
            # ACTAG,AT,.,.,.,GCC
            for idx, (off_target_char, on_target_char) in enumerate(zip(off_target_aln, on_target_aln)):
                if on_target_char == off_target_char:
                    if first_base_is_del_in_ontarget: # if we finished a run of starting deletions, add this base to the previous char
                        this_chars_for_plot[-1] += off_target_char
                        first_base_is_del_in_ontarget = False
                    else:
                        this_chars_for_plot.append(".")
                elif on_target_char == '-': # if gap in ontarget (insertion in this guide/offtarget), add this base to the previous base
                    if idx == 0:
                        first_base_is_del_in_ontarget = True
                        this_chars_for_plot.append(off_target_char)
                    else:
                        if this_chars_for_plot[-1] == '.':
                            this_chars_for_plot[-1] = on_target_aln[idx-1]
                        this_chars_for_plot[-1] += off_target_char
                else: # mismatch
                    if first_base_is_del_in_ontarget: # if we finished a run of starting deletions, add this base to the previous char
                        this_chars_for_plot[-1] += off_target_char
                        first_base_is_del_in_ontarget = False
                    else:
                        this_chars_for_plot.append(off_target_char)

            post_len = len(''.join(this_chars_for_plot))
            if post_len != len(off_target_aln):
                print('offtarget lengths do not match (' + str(post_len) + ' vs ' + str(len(off_target_aln)) + ')') 
            if len(this_chars_for_plot) != len(on_target_aln.replace('-', '' )):
                print('ontarget lengths do not match (' + str(len(this_chars_for_plot)) + ' vs ' + str(len(on_target_aln.replace('-', '')) + ')'))

            this_id = row['guide_id']
            names_for_plot.append(this_id)
            seqs_for_plot.append((this_chars_for_plot, row['guide_pam']))
            aln_ontarget_seqs.append(on_target_aln)
            aln_offtarget_seqs.append(off_target_aln)

        else: #this is the on-target
            last_guide_name = row['ontarget_name']
            last_guide_seq = row['ontarget_sequence']
            this_seq = row['guide_seq_no_gaps']
            this_id = row['guide_id']
            names_for_plot.append(this_id)
            seqs_for_plot.append((list(this_seq), row['guide_pam']))
            aln_ontarget_seqs.append(this_seq)
            aln_offtarget_seqs.append(this_seq)

    max_guide_len = 0
    max_pam_len = 0
    for guide_chars, pam_seq in seqs_for_plot:
        if len(guide_chars) > max_guide_len:
            max_guide_len = len(guide_chars)
        if len(pam_seq) > max_pam_len:
            max_pam_len = len(pam_seq)
        
    padded_seqs_for_plot = []
    for guide_chars, pam_seq in seqs_for_plot:
        padded_guide_chars = [' '] * (max_guide_len - len(guide_chars)) + guide_chars
        padded_pam_seq = pam_seq + ' ' * (max_pam_len - len(pam_seq))
        padded_seqs_for_plot.append(padded_guide_chars + list(padded_pam_seq))

    guide_df['padded_seqs_for_plot'] = padded_seqs_for_plot
    # these are for debugging - can be displayed as a table
    # guide_df['aln_ontarget_seqs'] = aln_ontarget_seqs
    # guide_df['aln_offtarget_seqs'] = aln_offtarget_seqs

    df_guides = pd.DataFrame(guide_df['padded_seqs_for_plot'].apply(lambda seq: [char for char in seq]).tolist())
    df_guides.index = guide_df['guide_id']
    return df_guides

def plot_guides_and_heatmap(guide_plot_df, df_data, col_to_plot, df_data_title, df_data_heat_max=None, df_data_heat_min=None, outfile_name=None, fig_height=24, fig_width=24, guide_names=None, seq_plot_ratio=1):
    """
    Plot a heatmap of guide sequences (left) and the data (right)
    The plot of guide sequences will contain 
        '-' if there is a base missing in the guide (relative to the on-target)
        '+' or 'AC' if there is a base added in the guide (relative to the on-target) 
    
    params:
    - guide_plot_df: a dataframe of guide sequences with each column showing a nucleotide position, and each row showing the bases for each guide
    - df_data: a dataframe of data to plot in a heatmap with rows corresponding to the rows in guide_plot_df
    - col_to_plot: the column suffix in df_data to plot (all columns with a name ending in this will be plotted)
    - df_data_title: the title of the heatmap
    - df_data_heat_max: the maximum value for the heatmap
    - df_data_heat_min: the minimum value for the heatmap
    - outfile_name: the name of the output file (if None, the plot will be displayed)
    - fig_height: the height of the figure
    - fig_width: the width of the figure
    - guide_names: a list of guide names to use for the y-axis of the heatmap (if None, the index of guide_plot_df will be used)
    - seq_plot_ratio: the ratio of the width of the guide sequence plot to the data plot (>1 means the seq plot is larger than the data plot)

    returns:
    - None
    """

    cols_to_plot = [col for col in df_data.columns if col_to_plot in col]
    if len(cols_to_plot) == 0:
        raise Exception('No columns found in df_data with suffix ' + col_to_plot + ' in columns: ' + str(list(df_data.columns)))
    df_to_plot = df_data[cols_to_plot]
    df_to_plot.columns = [col.replace('_' + col_to_plot, '') for col in cols_to_plot]

    # Create a custom color palette for the letters
    color_mapping = {'A': '#E3EFA9', 'T': '#CCCFE0', 'C': '#FBC6C6', 'G': '#FCE588', '.': '#F4F4F6', '-':'#BEC0C6', '+':'#BEC0C6', 'default': 'gray'}
    #color_mapping = {'A': '#90adc6', 'T': '#e9eaec', 'C': '#fad02c', 'G': '#333652', '.': '#c8df52', '-':'#F67E7D',  'default': 'gray'}
    #color_mapping = {'A': '#FFB7B2', 'T': '#C3CDE6', 'C': '#E2F0CB', 'G': '#B5EAD7', '.': '#FFDAC1', 'default': 'gray'}
    unique_letters = color_mapping.keys()
    colors = [color_mapping[letter] for letter in unique_letters]
    cmap = ListedColormap(colors)


    # Map letters to integers for heatmap
    letter_to_int = {letter: i for i, letter in enumerate(unique_letters)}
    def get_int_from_letter(letter):
        if letter in letter_to_int:
            return letter_to_int[letter]
        if len(letter) > 1:
            return letter_to_int['+']
        return letter_to_int['default']
    df_mapped = guide_plot_df.applymap(lambda x: get_int_from_letter(x))

    # Create the figure and gridspec
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(2, 2, height_ratios=[10, 1], width_ratios=[seq_plot_ratio, 1], hspace=0, wspace=0.05)

    # Create the first heatmap (guide_seqs and mismatches)
    ax1 = fig.add_subplot(gs[0, 0])
    sns.heatmap(df_mapped, annot=guide_plot_df, fmt='', cmap=cmap, cbar=False, ax=ax1, xticklabels=False, vmin=0, vmax=len(unique_letters)-1)
    #replace y tick lables with guide names
    if guide_names is not None:
        ax1.set_yticklabels(guide_names, rotation=0)
    ax1.set_title('Guide sequences')

    # Create the second heatmap (continuous data from df_data)
    ax2 = fig.add_subplot(gs[0, 1])
    data_annot = True
    if df_to_plot.shape[0] > 20:
        data_annot = False
    if df_data_heat_max is None:
        df_data_heat_max = df_to_plot.max(skipna=True).max(skipna=True)
    if df_data_heat_min is None:
        df_data_heat_min = df_to_plot.min(skipna=True).min(skipna=True)
    ax_heat = sns.heatmap(df_to_plot, annot=data_annot, cmap='Blues', ax=ax2, yticklabels=False, vmin=df_data_heat_min, vmax=df_data_heat_max)
    ax_heat.collections[0].cmap.set_bad('0.7')
    ax2.set_title(df_data_title)
    ax2.set_yticks([])

    # Create a custom legend under the first heatmap
    ax_legend = fig.add_subplot(gs[1, 0])
    ax_legend.axis('off')
    patches = [plt.plot([],[], marker="s", ms=10, ls="", mec=None, color=color_mapping[letter], 
                label="{:s}".format(letter) )[0]  for letter in color_mapping]
    ax_legend.legend(handles=patches, loc='center', title="Guide nucleotides", ncol=len(color_mapping), bbox_to_anchor=(0.5, 1))

    gs.tight_layout(fig)

    if outfile_name is not None:
        plt.savefig(outfile_name)
    else:
        plt.show()

    plt.close(fig)

def replot(reordered_guide_file, output_folder=None, file_prefix=None, name_column=None, fig_height=24, fig_width=24, seq_plot_ratio=1):
    """
    Replot a completed analysis using a reordered guide file

    params:
    - reordered_guide_file: the reordered guide file based on (having the same columns and layout as) a previously-completed aggregated_stats_all.txt file
    - output_folder: the output folder for the output plots
    - file_prefix: the prefix for the output plots
    - name_column: the column to use as the displayed name for each sample in the plot
    - fig_height: the height of the figure
    - fig_width: the width of the figure
    - seq_plot_ratio: the ratio of the width of the guide sequence plot to the data plot (>1 means the seq plot is larger than the data plot)
    """

    if output_folder is not None and file_prefix is not None:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder, exist_ok=True)
        output_root = output_folder + "/" + file_prefix + "."
    elif file_prefix is not None:
        output_root = file_prefix + "."
    elif output_folder is not None:
        output_root = output_folder + "/CRISPRessoSea."
    else:
        output_root = "CRISPRessoSea."
    reordered_guide_df = pd.read_csv(reordered_guide_file, sep="\t")
    guide_plot_df = create_guide_df_for_plotting(reordered_guide_df)

    if name_column is None:
        guide_plot_df.index = reordered_guide_df['guide_chr'] + ':' + reordered_guide_df['guide_pos'].astype(str) + ' ' + reordered_guide_df['guide_name']
    else:
        if name_column not in reordered_guide_df.columns:
            raise Exception('Name column ' + name_column + ' not found in reordered guide columns. Columns are: ' + str(reordered_guide_df.columns))
        guide_plot_df.index = reordered_guide_df[name_column]

    print('Plotting for ' + str(len(reordered_guide_df)) + ' guides')

    col_to_plot = 'highest_a_g_pct'
    col_title = 'Highest A/G %'
    plot_suffix = 'highest_a_g_pct.pdf'
    plot_guides_and_heatmap(guide_plot_df, reordered_guide_df, col_to_plot, col_title, outfile_name=output_root + plot_suffix, fig_height=fig_height, fig_width=fig_width, seq_plot_ratio=seq_plot_ratio)

    col_to_plot = 'highest_c_t_pct'
    col_title = 'Highest C/T %'
    plot_suffix = 'highest_c_t_pct.pdf'
    plot_guides_and_heatmap(guide_plot_df, reordered_guide_df, col_to_plot, col_title, outfile_name=output_root + plot_suffix, fig_height=fig_height, fig_width=fig_width, seq_plot_ratio=seq_plot_ratio)

    col_to_plot = 'highest_indel_pct'
    col_title = 'Highest Indel %'
    plot_suffix = 'highest_indel_pct.pdf'
    plot_guides_and_heatmap(guide_plot_df, reordered_guide_df, col_to_plot, col_title, outfile_name=output_root + plot_suffix, fig_height=fig_height, fig_width=fig_width, seq_plot_ratio=seq_plot_ratio)


# main entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process multiple pooled sequencing runs')
    subparsers = parser.add_subparsers(dest='subcommand', help='Process a new run or Replot a completed run')
    process_parser = subparsers.add_parser('Process', help='Process a new set of samples', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    process_parser.add_argument('-o', '--output_folder', help='Output folder', default=None)
    process_parser.add_argument('-g','--guide_file', help='Guide file - list of guides with one guide per line', required=True)
    process_parser.add_argument('-s', '--sample_file', help='Sample file - list of samples with one sample per line', required=True)
    process_parser.add_argument('-x', '--genome_file', help='Bowtie2-indexed genome file - files ending in .bt2 must be present in the folder.', required=True)
    process_parser.add_argument('-p', '--n_processes', help='Number of processes to use', type=int, default=8)
    process_parser.add_argument('--allow_unplaced_chrs', help='Allow regions on unplaced chromosomes (chrUn, random, etc). By default, regions on these chromosomes are excluded. If set, regions on these chromosomes will be included.', action='store_true')
    process_parser.add_argument('--plot_only_complete_guides', help='Plot only guides with all values. If not set, all guides will be plotted.', action='store_true')
    process_parser.add_argument('--min_amplicon_coverage', help='Minimum number of reads to cover a location for it to be plotted. Otherwise, it will be set as NA', default=10, type=int)
    process_parser.add_argument('--sort_based_on_mismatch', help='Sort guides based on mismatch count. If true, the on-target will always be first', action='store_true')
    process_parser.add_argument('--allow_guide_match_to_other_region_loc', help='If true, guides can match to regions even if the guide chr:start is not in that region (e.g. if the guide sequence is found in that region). If false/unset, guides can only match to regions matching the guide chr:start position. This flag should be set if the genome for guide design was not the same as the analysis genome.', action='store_true')

    plot_parser = subparsers.add_parser('Replot', help='Replot completed analysis using reordered sample table', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    plot_parser.add_argument('-o', '--output_folder', help='Output folder', default=None)
    plot_parser.add_argument('-p', '--file_prefix', help='File prefix for output files', default='CRISPRessoSea')
    plot_parser.add_argument('-f','--reordered_guide_file', help='Reordered guide file - made by reordering rows from aggregated_stats_all.txt', required=True)
    plot_parser.add_argument('-n','--name_column', help='Column name to set as the displayed name for each sample in the plot', default=None)
    plot_parser.add_argument('--fig_width', help='Width of the figure', default=24, type=int)
    plot_parser.add_argument('--fig_height', help='Height of the figure', default=24, type=int)
    plot_parser.add_argument('--seq_plot_ratio', help='Ratio of the width of the sequence plot to the data plot (>1 means the seq plot is larger than the data plot)', default=1, type=float)

    args = parser.parse_args()
    
    if args.subcommand == 'Replot':
        if not os.path.isfile(args.reordered_guide_file):
            raise Exception('Reordered guide file is not found at "' + args.reordered_guide_file + '". Use an aggregated_stats_all.txt file from a previous run as a template.')
        replot(args.reordered_guide_file, args.output_folder, args.file_prefix, name_column=args.name_column, fig_height=args.fig_height, fig_width=args.fig_width, seq_plot_ratio=args.seq_plot_ratio)

    elif args.subcommand == 'Process':
        output_folder = args.output_folder
        if output_folder is None:
            output_folder = 'CRISPRessoSea_output_on_' + os.path.basename(args.sample_file)
        if not output_folder.endswith('/'):
            output_folder += '/'
        os.makedirs(output_folder, exist_ok=True)

        if not os.path.isfile(args.sample_file):
            raise Exception('Sample file not found at ' + args.sample_file)

        if not os.path.isfile(args.guide_file):
            raise Exception('Guide file not found at ' + args.guide_file)

        process_pools(sample_file=args.sample_file, guide_file=args.guide_file, genome_file=args.genome_file,
            output_folder=output_folder, n_processors=args.n_processes, allow_unplaced_chrs=args.allow_unplaced_chrs, 
            plot_only_complete_guides=args.plot_only_complete_guides, min_amplicon_coverage=args.min_amplicon_coverage,
            sort_based_on_mismatch=args.sort_based_on_mismatch, allow_guide_match_to_other_region_loc=args.allow_guide_match_to_other_region_loc)
    else:
        raise Exception('Please run with subcommand "Process" (to process initial analysis) or "Replot" (to replot finished analysis)')
    print('Finished')
