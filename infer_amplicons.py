import argparse
import subprocess
import pandas as pd
import os
import sys
from collections import defaultdict

from CRISPResso2 import CRISPRessoShared

experiment_name = 'test1'
guide_file = '02_parse_from_idt/ordered_guides.txt.forCRISPResso.txt'
fastq_r1 = 'test1_R1.fastq'
fastq_r2 = 'test1_R2.fastq'
genome_file = 'genome.fasta'

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))

def merge_locations(crispresso_pooled_genome_folder, skip_bad_chrs=True, min_amplicon_len=50, debug=False):
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
        
    for index, region_row in region_df_sub.iterrows():
        region_chr = region_row.chr_id
        region_start = region_row.start
        region_end = region_row.end
        region_count = region_row['number of reads']
        region_name = "_".join([region_chr, str(region_start), str(region_end)])

        region_info[region_name] = {'chr': region_chr, 'start': region_start, 'end': region_end, 'region_count': region_count}
        
        if skip_bad_chrs:
            if 'chrUn' in region_chr or 'random' in region_chr:
                continue
        region_len = region_end - region_start
        if region_len < min_amplicon_len:
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


def run_demux(experiment_name, guide_file, fastq_r1, fastq_r2, genome_file, guides_have_pam=False, n_processors=8, skip_bad_chrs=True):
    # Run CRISPResso

    if fastq_r1 == fastq_r2:
        raise Exception('R1 and R2 files must be different')
    
    outfile_root = experiment_name
    CRISPResso_output_folder = "CRISPRessoPooled_on_"+experiment_name
    command = 'CRISPRessoPooled -x ' + genome_file + ' -r1 ' + fastq_r1 + ' -r2 ' + fastq_r2 + ' -n ' + experiment_name + ' -p ' + str(n_processors) + ' --no_rerun --keep_intermediate --suppress_plots'

    crispresso_run_is_complete = False
    if os.path.exists(CRISPResso_output_folder):
        crispresso_pooled_info = CRISPRessoShared.load_crispresso_info(crispresso_info_file_path=CRISPResso_output_folder + '/CRISPResso2Pooled_info.json')
        if 'demultiplexing_genome_only_regions' in crispresso_pooled_info['running_info']['finished_steps']:
            crispresso_run_is_complete = True
        print('CRISPResso output folder already exists, skipping CRISPResso run')
    if not crispresso_run_is_complete:
        print('Aligning reads to the genome to find aligned locations')
        print('Running command ' + str(command))
        subprocess.run(command, shell=True, check=True)

    guide_seq_ids = []
    guide_sequences = []
    unique_guide_seq_ids = []
    guide_annotations = {}
    with open(guide_file) as f:
        read_guide_seq_counts = 1
        for line in f:
            line_els = line.strip().split("\t")
            guide_seq = line_els[0]
            this_guide_seq_id = str(read_guide_seq_counts) + "_" + guide_seq
            this_guide_anno = None
            if len(line_els) > 1:
                this_guide_anno = line_els[1]
            if len(guide_seq) > 0:
                guide_seq_ids.append(this_guide_seq_id)
                if guide_seq not in guide_sequences: # guide is unique (haven't seen it before)
                    unique_guide_seq_ids.append(this_guide_seq_id)
                guide_annotations[this_guide_seq_id] = (guide_seq, this_guide_anno)
            read_guide_seq_counts += 1

    print('got ' + str(len(unique_guide_seq_ids)) + '/' + str(len(guide_seq_ids)) + ' unique guide sequences')

    merged_regions, merged_regions_infos = merge_locations(CRISPResso_output_folder, skip_bad_chrs=skip_bad_chrs)
    for region in merged_regions:
        region_chr = merged_regions_infos[region]['chr']
        region_start = merged_regions_infos[region]['start']
        region_end = merged_regions_infos[region]['end']

        region_seq_output = subprocess.check_output(
            '%s faidx %s %s:%d-%d'%('samtools',genome_file + ".fa",region_chr,region_start,region_end-1),shell=True).decode(sys.stdout.encoding)
        region_seq = "".join(region_seq_output.split("\n")[1:])
        merged_regions_infos[region]['seq'] = region_seq


    match_count = 0
    unmatch_count = 0

    target_matches = defaultdict(list)
    all_region_matches = defaultdict(list)
    for guide_seq_id in guide_seq_ids:
        guide_is_matched = False
        guide_seq, guide_anno = guide_annotations[guide_seq_id]
        for region in merged_regions:
            region_info = merged_regions_infos[region]
            region_chr = region_info['chr']
            region_start = region_info['start']
            region_end = region_info['end']
            region_seq = region_info['seq']
            
            is_match = False
            target_seq_match = False

            if guide_seq.upper() in region_seq.upper():
                is_match = True
                guide_is_matched = True
            elif reverse_complement(guide_seq).upper() in region_seq.upper():
                is_match = True
                guide_is_matched = True

            if is_match:
                target_matches[guide_seq_id].append(region)
                all_region_matches[region].append(guide_seq_id)
                target_is_matched = True

        if not guide_is_matched:
            print('no match for guide: ' + guide_seq_id + ": " + str(guide_seq))
            unmatch_count += 1
        else:
            match_count += 1

    out_file = outfile_root + '.matches.txt'
    region_matches = {} # region_name -> guide_seq
    printed_regions = {}
    with open(out_file, 'w') as fout:
        fout.write('guide_seq_id\tguide_seq\tguide_seq_annot\tchr\tstart\tend\tread_count\tregion_seq\n')
        for guide_seq_id in guide_seq_ids:
            guide_seq, guide_anno = guide_annotations[guide_seq_id]
            if guide_seq_id in target_matches:
                matched_region_count = len(target_matches[guide_seq_id])
                region_name = target_matches[guide_seq_id][0]
                region_info = merged_regions_infos[region_name]
                for other_match in target_matches[guide_seq_id]:
                    other_match_info = merged_regions_infos[other_match]
                    if other_match_info['region_count'] > region_info['region_count']:
                        region_name = other_match
                        region_info = other_match_info

                if region_name in region_matches:
                    print('Warning: multiple guides matched to the same region ' + region_name + '(' + guide_seq_id + ' and ' + region_matches[region_name] + ')')
                region_matches[region_name] = guide_seq_id
                region_chr = region_info['chr']
                region_start = region_info['start']
                region_end = region_info['end']
                region_seq = region_info['seq']

                fout.write('\t'.join([str(x) for x in [guide_seq_id, guide_seq, guide_anno, region_chr, region_start, region_end, region_info['region_count'], region_seq]]) + '\n')
                region_name = "_".join([region_chr, str(region_start), str(region_end)])
                crispresso_guide = guide_seq
                if guides_have_pam:
                    crispresso_guide = guide_seq[:-3]

                this_region_name = region_name
                if guide_annotations[guide_seq_id] is not None:
                    this_region_name = region_name + '_' + guide_anno
                if region_seq not in printed_regions:
                    printed_regions[region_seq] = region_name
                else:
                    print('For guide ' + guide_seq + ' - skipping adding a CRISPResso region ' + region_name + ' because its sequence is a duplicate of ' + printed_regions[region_seq])

    with open(outfile_root + '.guide_matches.txt', 'w') as fout:
        fout.write('guide_seq_id\tguide_seq\tguide_anno\tmatched_region_count\tmatched_region_names\tchr\tstart\tend\tregion_count\tregion_seq\n')
        for guide_seq_id in guide_seq_ids:
            guide_seq, guide_anno = guide_annotations[guide_seq_id]
            if guide_seq_id in target_matches:
                matched_region_count = len(target_matches[guide_seq_id])
                region_name = target_matches[guide_seq_id][0]
                region_info = merged_regions_infos[region_name]
                region_chr = region_info['chr']
                region_start = region_info['start']
                region_end = region_info['end']
                region_seq = region_info['seq']
                fout.write('\t'.join([str(x) for x in [guide_seq_id, guide_seq, guide_anno, matched_region_count, ",".join(target_matches[guide_seq_id]), region_chr, region_start, region_end, region_info['region_count'], region_seq]]) + '\n')
            else:
                fout.write('\t'.join([str(x) for x in [guide_seq_id, guide_seq, guide_anno, 0, '', '', '', '', '', '']]) + '\n')


    crispresso_output_file = outfile_root + '.CRISPRessoPooledRegions.txt'

    all_region_output_file = outfile_root + '.all_regions.txt'
    with open(all_region_output_file,'w') as rout, open(crispresso_output_file,'w') as cout:
        rout.write('type\tchr\tstart\tend\tread_count\tseq\ttarget\tguide_seq\tguide_anno\n')
        for region in merged_regions:
            region_info = merged_regions_infos[region]
            region_chr = region_info['chr']
            region_start = region_info['start']
            region_end = region_info['end']
            region_seq = region_info['seq']
            region_name = "_".join([region_chr, str(region_start), str(region_end)])
            region_target = 'NA'
            region_anno = 'NA'
            guide_seq = 'NA'
            guide_anno = 'NA'
            if region in all_region_matches:
                this_guide_ids = all_region_matches[region]
                if len(this_guide_ids) > 1:
                    print('Warning: multiple guides matched to the same region ' + region + ': ' + str(this_guide_ids))
                this_guide_id = this_guide_ids[0]
                guide_seq, guide_anno = guide_annotations[this_guide_id]
                this_region_name = region_name
                if guide_annotations[guide_seq_id] is not None:
                    this_region_name = region_name + '_' + guide_anno
                cout.write("\t".join([this_region_name,region_seq,guide_seq]) + '\n')
            rout.write('\t'.join(str(x) for x in ['region', region_chr, region_start, region_end, region_info['region_count'], region_seq, region_target, guide_seq, guide_anno]) + '\n')
        
            
    print('Matched ' + str(match_count) + "/" + str(len(unique_guide_seq_ids)) + ' targets from metadata to regions from sequencing analysis. Wrote matches to ' + out_file)

    # Run CRISPResso again on the subset, this time including guide sequences
    command = 'CRISPRessoPooled -f ' + crispresso_output_file + ' -x ' + genome_file + ' -n ' + experiment_name + '_polished --min_reads_to_use_region 1 --default_min_aln_score 20 ' + \
                ' --base_editor_output --quantification_window_center -10 --quantification_window_size 16' + \
                ' -p ' + str(n_processors) + ' --no_rerun  --exclude_bp_from_left 0 --exclude_bp_from_right 0 --plot_window_size 10 --aligned_pooled_bam ' + CRISPResso_output_folder + '/' + experiment_name + '_GENOME_ALIGNED.bam' 
    print('running command ' + str(command))
    subprocess.run(command, shell=True, check=True)

# main entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Infer amplicon sequences from unannotated read file')
    parser.add_argument('-n', '--experiment_name', help='Name of the experiment', required=True)
    parser.add_argument('-g','--guide_file', help='Guide file - list of guides with one guide per line', required=True)
    parser.add_argument('-r1', '--fastq_r1', help='Read 1 fastq file', required=True)
    parser.add_argument('-r2', '--fastq_r2', help='Read 2 fastq file', required=True)
    parser.add_argument('-x', '--genome_file', help='Bowtie2-indexed Genome file', required=True)
    parser.add_argument('-p', '--n_processes', help='Number of processes to use', type=int, default=8)
    parser.add_argument('--guides_have_pam', help='Provided guides have PAM sequences (if unset, guides are provided without PAM sequences)', action='store_true')

    args = parser.parse_args()

    run_demux(args.experiment_name, args.guide_file, args.fastq_r1, args.fastq_r2, args.genome_file, args.guides_have_pam)
    print('Finished')
