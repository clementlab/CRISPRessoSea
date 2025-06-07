import argparse
import glob
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import subprocess
import sys
import scipy.stats as sp
import statsmodels.api as sm


from CRISPResso2 import CRISPResso2Align, CRISPRessoMultiProcessing, CRISPRessoShared
from CRISPResso2.CRISPRessoReports import CRISPRessoReport

from collections import defaultdict
from copy import deepcopy
from enum import Enum
from jinja2 import Environment, FileSystemLoader, make_logging_undefined
from matplotlib.colors import ListedColormap
from statsmodels.api import GLM, families
from statsmodels.stats.multitest import multipletests
from statsmodels.discrete.discrete_model import NegativeBinomial

mpl.rcParams["pdf.fonttype"] = 42

__version__ = "0.1.2"

C2PRO_INSTALLED = False

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(CRISPRessoShared.LogStreamHandler())

error = logger.critical
warn = logger.warning
debug = logger.debug
info = logger.info


def process_pools(
    args,
    sample_file,
    target_file,
    genome_file,
    output_folder,
    gene_annotations=None,
    n_processors=8,
    crispresso_quantification_window_center=-3,
    crispresso_quantification_window_size=1,
    crispresso_base_editor_output=False,
    crispresso_default_min_aln_score=20,
    crispresso_plot_window_size=20,
    allow_unplaced_chrs=False,
    plot_only_complete_targets=False,
    min_amplicon_coverage=100,
    sort_based_on_mismatch=False,
    allow_target_match_to_other_region_loc=False,
    top_percent_cutoff=0.2,
    min_amplicon_len=50,
    fail_on_pooled_fail=False,
    plot_group_order=None,
    sig_method_parameters=None,
):
    """
    Discover amplicons and analyze sample editing at amplicons. This is the main entrypoint for this program.

    params:
    - sample_file: path to the sample file with headers: Name, group, fastq_r1, fastq_r2 (group is always optional, fastq_r2 is optional for single-end reads)
    - target_file: path to the target file with headers: Guide, Sequence, PAM, #MM, Locus (may include other columns)
    - genome_file: path to the genome file (not including the trailing .fa)
    - output_folder: path to the output folder to produce data
    - gene_annotations: the path to the gene annotations file. This function expects a column for:
        chromosome ('chrom' or 'chr'),
        start ('txStart','start', or 'Start'),
        end ('txEnd', 'end', or 'End'), and
        gene name ('name' and/or 'name2')
    - n_processors: number of processors to use for multiprocessing
    - crispresso_quantification_window_center: the window to use for CRISPResso2 quantification. This will be used for CRISPResso2 subruns and corresponds to the parameter --quantification_window_center. If None, the default value of -3 will be used.
    - crispresso_quantification_window_size: the window size to use for CRISPResso2 quantification. This will be used for CRISPResso2 subruns and corresponds to the parameter --quantification_window_size. If None, the default value of 1 will be used.
    - crispresso_base_editor_output: whether to output base editor plots. If True, base editor plots will be output. If False, base editor plots will not be output.
    - crispresso_default_min_aln_score: the minimum alignment score to use for CRISPResso2. This will be used for CRISPResso2 subruns and corresponds to the parameter --default_min_aln_score. If None, the default value of 20 will be used.
    - crispresso_plot_window_size: the window size to use for CRISPResso2 plots. This will be used for CRISPResso2 subruns and corresponds to the parameter --plot_window_size. If None, the default value of 20 will be used.
    - allow_unplaced_chrs: whether to allow regions to be identified on unplaced chromosomes (chrUn, random, etc).
    - plot_only_complete_targets: whether to plot only targets with data in all samples
    - min_amplicon_coverage: the minimum number of reads required to consider an amplicon
    - sort_based_on_mismatch: if true, targets are sorted based on mismatch count. If false, targets are presented in file order
    - allow_target_match_to_other_region_loc: if true, targets can match regions based on sequence even if not in the same pos.
    - top_percent_cutoff: the top percent of aligned regions (by region read depth) to consider in finding non-overlapping regions. This is a float between 0 and 1. For example, if set to 0.2, the top 20% of regions (by read depth) will be considered.
    - min_amplicon_len: the minimum length of an amplicon to consider in finding non-overlapping regions. Amplicons shorter than this will be ignored.
    - fail_on_pooled_fail: if true, fail if any pooled CRISPResso run fails. Otherwise, continue even if sub-CRISPResso commands fail.
    - plot_group_order: the order of groups to plot. If None, the default order will be used. This is a list of strings.
    - sig_method_parameters: the parameters to use for the significance method. 
    """

    log_filename = os.path.join(output_folder, "CRISPRessoSea_RUNNING_LOG.txt")
    logger.addHandler(logging.FileHandler(log_filename))
    status_handler = CRISPRessoShared.StatusHandler(
        os.path.join(output_folder, "CRISPRessoSea_status.json")
    )
    logger.addHandler(status_handler)

    with open(log_filename, "w+") as outfile:
        outfile.write("[Command used]:\n%s\n\n[Execution log]:\n" % " ".join(sys.argv))

    crispressoSea_info_file = os.path.join(output_folder, "CRISPRessoSea_info.json")

    crispresso2_info = {
        "running_info": {
            "version": __version__,
            "log_filename": os.path.basename(log_filename),
            "args": deepcopy(args),
        },
        "results": {
            "summary_stats": {},
            "alignment_stats": {},
            "completed_sea_arr": [],
            "failed_sea_arr": [],
            "failed_sea_arr_desc": {},
            "sea_input_names": {},
            "general_plots": {
                "summary_plot_names": [],
                "summary_plot_paths": {},
                "summary_plot_titles": {},
                "summary_plot_labels": {},
                "summary_plot_datas": {},
                "summary_plot_divs": {},
            },
            "paths": {"graphs": []},
            "samples": [],
        },
    }

    CRISPRessoShared.write_crispresso_info(
        crispressoSea_info_file,
        crispresso2_info,
    )

    # Load sample file
    sample_df = parse_sample_file(sample_file)

    # pull out first sample name which we will use to find amplicon locations
    first_sample = sample_df.iloc[0, :]
    first_sample_name = first_sample["Name"]
    first_sample_r1 = first_sample["fastq_r1"]
    first_sample_r2 = None
    if "fastq_r2" in first_sample:
        first_sample_r2 = first_sample["fastq_r2"]
    if genome_file.endswith(".fa"):
        genome_file = genome_file[:-3]

    # Check fastq files for existence
    for idx, row in sample_df.iterrows():
        sample_name = row["Name"]
        sample_r1 = row["fastq_r1"]
        sample_r2 = None
        if "fastq_r2" in row and row["fastq_r2"] is not None:
            sample_r2 = row["fastq_r2"]
        if not os.path.exists(sample_r1):
            raise Exception(
                "Fastq R1 file for sample " + sample_name + " not found at " + sample_r1
            )
        if sample_r2 is not None and not os.path.exists(sample_r2):
            raise Exception(
                "Fastq R2 file for sample " + sample_name + " not found at " + sample_r2
            )

    summary_output_folder = os.path.join(output_folder, "summary/")
    if not os.path.exists(summary_output_folder):
        os.makedirs(summary_output_folder, exist_ok=True)

    crispresso_output_folder = os.path.join(output_folder, "CRISPResso_output/")
    if not os.path.exists(crispresso_output_folder):
        os.makedirs(crispresso_output_folder, exist_ok=True)

    # Define variables referenced later but not defined previously
    completed_samples = []
    failed_samples = {}
    display_names = {}

    merged_region_info_file = run_initial_demux(
        first_sample_name,
        first_sample_r1,
        first_sample_r2,
        genome_file,
        output_folder=crispresso_output_folder,
        n_processors=n_processors,
        crispresso_quantification_window_center=crispresso_quantification_window_center,
        crispresso_quantification_window_size=crispresso_quantification_window_size,
        crispresso_base_editor_output=crispresso_base_editor_output,
        crispresso_default_min_aln_score=crispresso_default_min_aln_score,
        crispresso_plot_window_size=crispresso_plot_window_size,
        allow_unplaced_chrs=allow_unplaced_chrs,
        top_percent_cutoff=top_percent_cutoff,
        min_amplicon_len=min_amplicon_len,
        fail_on_pooled_fail=fail_on_pooled_fail,
    )

    crispresso_region_file, target_df, region_df = make_target_region_assignments(
        merged_region_info_file,
        target_file,
        output_folder,
        sort_based_on_mismatch,
        allow_target_match_to_other_region_loc=allow_target_match_to_other_region_loc,
    )
    if gene_annotations is not None:
        target_df = add_region_annotations_to_target_df(target_df, gene_annotations)
    else:
        target_df["region_anno"] = ""

    aggregated_stats = target_df.copy()
    aggregated_stats.set_index("target_id", inplace=True)
    sample_df["CRISPRessoPooled_output_folder"] = "NA"

    start_pct_complete = 40
    end_pct_complete = 90
    current_pct_complete = start_pct_complete
    step_pct_complete = (end_pct_complete - start_pct_complete) / len(sample_df)

    for sample_idx, sample_row in sample_df.iterrows():
        sample_name = sample_row["Name"]
        current_pct_complete += step_pct_complete
        info('Processing sample ' + sample_name, {'percent_complete': current_pct_complete})
        sample_r1 = sample_row["fastq_r1"]
        sample_r2 = None
        if "fastq_r2" in sample_row:
            sample_r2 = sample_row["fastq_r2"]
        this_pooled_run = run_crispresso_with_assigned_regions(
            sample_name,
            sample_r1,
            sample_r2,
            genome_file,
            crispresso_region_file,
            crispresso_output_folder,
            n_processors=n_processors,
            crispresso_quantification_window_center=crispresso_quantification_window_center,
            crispresso_quantification_window_size=crispresso_quantification_window_size,
            crispresso_base_editor_output=crispresso_base_editor_output,
            crispresso_default_min_aln_score=crispresso_default_min_aln_score,
            crispresso_plot_window_size=crispresso_plot_window_size,
            fail_on_pooled_fail=fail_on_pooled_fail,
        )
        sample_df.loc[sample_idx, "CRISPRessoPooled_output_folder"] = this_pooled_run
        if this_pooled_run is not None:
            target_summary_file, target_summary_df = analyze_run(
                os.path.join(summary_output_folder, sample_name),
                this_pooled_run,
                target_df,
                region_df,
            )
            sample_df.loc[sample_idx, "target_summary_file"] = target_summary_file

            target_summary_df.columns = ["target_id", "target_label"] + [
                sample_name + "_" + x for x in target_summary_df.columns[2:]
            ]
            target_summary_df.drop(columns=["target_label"], inplace=True)

            aggregated_stats = pd.merge(aggregated_stats, target_summary_df, how="left", on="target_id")

            # for the columns in target_summary_df, if tot_reads is less than min_amplicon_coverage, set all other columns to NA
            for col in target_summary_df.columns:
                if col.startswith(sample_name) and col != sample_name + "_tot_reads":
                    aggregated_stats.loc[aggregated_stats[sample_name + "_tot_reads"] < min_amplicon_coverage, col] = np.nan
            completed_samples.append(sample_name)

    aggregated_stats.sort_values(by="sort_index", inplace=True)
    aggregated_stats.fillna('NA').to_csv(os.path.join(output_folder, "aggregated_stats_all.txt"), sep="\t", index=False)

    aggregated_stats["target_id"] = aggregated_stats["target_id"].apply(
        lambda x: x.split(" ")[-1]
    )

    if plot_only_complete_targets:
        aggregated_stats_good = aggregated_stats.dropna()
        info(
            "Plotting for "
            + str(len(aggregated_stats_good))
            + "/"
            + str(len(aggregated_stats))
            + " targets after removing targets missing data in any samples"
        )
    else:
        aggregated_stats_good = aggregated_stats
        info(
            "Plotting for "
            + str(len(aggregated_stats_good))
            + " targets"
        )

    target_plot_df = create_target_df_for_plotting(aggregated_stats_good)
    target_plot_df.index = aggregated_stats_good["target_name"] # this column isn't in the target_plot_df

    crispresso2_info = create_plots(
        data_df=aggregated_stats_good,
        sample_df=sample_df,
        target_plot_df=target_plot_df,
        output_folder=output_folder,
        file_prefix=None,
        crispresso2_info=crispresso2_info,
        plot_group_order=plot_group_order,
        sig_method_parameters=sig_method_parameters
    )
    # Update crispresso2_info fields that depend on sample_df
    crispresso2_info["results"]["samples"] = [
        (row["Name"], row.get("group", "")) for _, row in sample_df.iterrows()
    ]
    crispresso2_info["results"]["completed_sea_arr"] = completed_samples
    crispresso2_info["results"]["failed_sea_arr"] = list(failed_samples.keys())
    crispresso2_info["results"]["failed_sea_arr_desc"] = failed_samples
    crispresso2_info["results"]["sea_input_names"] = display_names
    crispresso2_info["results"]["summary_stats"] = aggregated_stats
    crispresso2_info["results"]["sea_input_names"] = {
        sample_name: sample_name for sample_name in sample_df["Name"]
    }


    crispresso2_info["results"]["sea_input_groups"] = {}
    for idx, row in sample_df.iterrows():
        sample_name = row["Name"]
        group = row.get("group", "")
        if group != '':
            crispresso2_info["results"]["sea_input_groups"][sample_name] = group

    _root = os.path.abspath(os.path.dirname(__file__))
    crispressoSea_report_file = os.path.join(
        output_folder, "output_crispresso_sea.html"
    )
    sea_folder = output_folder

    make_sea_report_from_folder(
        crispressoSea_report_file, crispresso2_info, sea_folder, _root, logger
    )

    CRISPRessoShared.write_crispresso_info(
            crispressoSea_info_file,
            crispresso2_info,
        )

    info('Analysis Complete!', {'percent_complete': 100})




def get_jinja_loader(root, logger):
    """
    Get the Jinja2 environment for rendering templates.
    """
    undefined_logger = make_logging_undefined(logger=logger)
    return Environment(
        loader=FileSystemLoader(os.path.join(root,'templates')),
        undefined=undefined_logger,
    )


# shim from https://github.com/pinellolab/CRISPResso2/blob/0232b0625d1867d268325a5f4162fb482ab6e0e4/CRISPResso2/CRISPRessoReports/CRISPRessoReport.py#L344C5-L344C22
def make_multi_report(
    run_names,
    failed_runs,
    failed_runs_desc,
    sub_html_files,
    crispresso_multi_report_file,
    crispresso_folder,
    _root,
    report_name,
    crispresso_tool,
    logger,
    window_nuc_pct_quilts=None,
    nuc_pct_quilts=None,
    window_nuc_conv_plots=None,
    nuc_conv_plots=None,
    summary_plots=None,
    compact_plots_to_show=None,
    allele_modification_heatmap_plot=None,
    allele_modification_line_plot=None,
):
    """
    Makes an HTML report for a run containing multiple crispresso runs

    Parameters:
    run_names (arr of strings): names of runs
    sub_html_files (dict): dict of run_name->file_loc
    crispresso_multi_report_file (string): path of file to write to
    report_name (string): description of report type to be shown at top of report
    crispresso_folder (string): absolute path to the crispresso output
    _root (string): absolute path to the crispresso executable
    summary_plots (dict): a dict with the following keys:
        names (list): list of plot names - keys for following dicts
        titles (dict): dict of plot_name->plot_title
        labels (dict): dict of plot_name->plot_label
        datas (dict): dict of plot_name->[(datafile_description, data_filename), ...]
    compact_plots_to_show (dict): name=>{'href': path to target(report) when user clicks on image, 'img': path to png image to show}
    allele_modification_heatmap_plot (dict): a dict with the following keys:
        names (list): list of plot names for heatmaps, keys for dicts below
        htmls (dict): dict of plot_name->HTML for the plot
        titles (dict): dict of plot_name->plot_title
        labels (dict): dict of plot_name->plot_label
        datas (dict): dict of plot_name->[(datafile_description, data_filename), ...]
    """

    def dirname(path):
        return os.path.basename(os.path.dirname(path))

    def fill_default(dictionary, key, default_type=list):
        if key not in dictionary:
            dictionary[key] = default_type()

    j2_env = get_jinja_loader(_root, logger)

    j2_env.filters["dirname"] = dirname
    if crispresso_tool == "sea":
        template = "seaReport.html"
    else:
        raise Exception('Cannot create report for tool "' + crispresso_tool + '"')

    crispresso_data_path = os.path.relpath(
        crispresso_folder,
        os.path.dirname(crispresso_multi_report_file),
    )
    if crispresso_data_path == ".":
        crispresso_data_path = ""
    else:
        crispresso_data_path += "/"

    if allele_modification_heatmap_plot is None:
        allele_modification_heatmap_plot = {}
    if allele_modification_line_plot is None:
        allele_modification_line_plot = {}
    dictionaries = [
        allele_modification_heatmap_plot,
        allele_modification_line_plot,
    ]
    keys_and_default_types = [
        ("names", list),
        ("htmls", dict),
        ("titles", list),
        ("labels", dict),
        ("datas", dict),
        ("divs", dict),
    ]
    for dictionary in dictionaries:
        for key, default_type in keys_and_default_types:
            fill_default(
                dictionary,
                key,
                default_type,
            )
    if summary_plots is None:
        summary_plots = {
            "names": [],
            "titles": [],
            "labels": [],
            "datas": [],
            "htmls": [],
        }

    for html in sub_html_files:
        sub_html_files[html] = crispresso_data_path + sub_html_files[html]
    with open(crispresso_multi_report_file, "w", encoding="utf-8") as outfile:
        outfile.write(
            CRISPRessoReport.render_template(
                template,
                j2_env,
                window_nuc_pct_quilts=(
                    [] if window_nuc_pct_quilts is None else window_nuc_pct_quilts
                ),
                nuc_pct_quilts=[] if nuc_pct_quilts is None else nuc_pct_quilts,
                window_nuc_conv_plots=(
                    [] if window_nuc_conv_plots is None else window_nuc_conv_plots
                ),
                nuc_conv_plots=[] if nuc_conv_plots is None else nuc_conv_plots,
                crispresso_data_path=crispresso_data_path,
                report_data={
                    "names": summary_plots["names"],
                    "titles": summary_plots["titles"],
                    "labels": summary_plots["labels"],
                    "datas": summary_plots["datas"],
                    "htmls": summary_plots["htmls"] if "htmls" in summary_plots else [],
                    "crispresso_data_path": crispresso_data_path,
                },
                run_names=run_names,
                failed_runs=failed_runs,
                failed_runs_desc=failed_runs_desc,
                sub_html_files=sub_html_files,
                report_name=report_name,
                compact_plots_to_show=(
                    [] if compact_plots_to_show is None else compact_plots_to_show
                ),
                allele_modification_heatmap_plot_names=allele_modification_heatmap_plot[
                    "names"
                ],
                allele_modification_heatmap_plot_htmls=allele_modification_heatmap_plot[
                    "htmls"
                ],
                allele_modification_heatmap_plot_titles=allele_modification_heatmap_plot[
                    "titles"
                ],
                allele_modification_heatmap_plot_labels=allele_modification_heatmap_plot[
                    "labels"
                ],
                allele_modification_heatmap_plot_datas=allele_modification_heatmap_plot[
                    "datas"
                ],
                allele_modification_heatmap_plot_divs=allele_modification_heatmap_plot[
                    "divs"
                ],
                allele_modification_line_plot_names=allele_modification_line_plot[
                    "names"
                ],
                allele_modification_line_plot_htmls=allele_modification_line_plot[
                    "htmls"
                ],
                allele_modification_line_plot_titles=allele_modification_line_plot[
                    "titles"
                ],
                allele_modification_line_plot_labels=allele_modification_line_plot[
                    "labels"
                ],
                allele_modification_line_plot_datas=allele_modification_line_plot[
                    "datas"
                ],
                allele_modification_line_plot_divs=allele_modification_line_plot[
                    "divs"
                ],
                C2PRO_INSTALLED=C2PRO_INSTALLED,
            )
        )


def make_sea_report_from_folder(
    crispressoSea_report_file, crispresso2_info, sea_folder, _root, logger
):
    """
    Makes a report for a CRIPSRessoSea run
    """
    sea_names = crispresso2_info["results"]["completed_sea_arr"]
    failed_runs = crispresso2_info["results"]["failed_sea_arr"]
    failed_runs_desc = crispresso2_info["results"]["failed_sea_arr_desc"]
    display_names = crispresso2_info["results"]["sea_input_names"]
    groups = crispresso2_info["results"]["sea_input_groups"]

    summary_plot_names = []
    if "summary_plot_names" in crispresso2_info["results"]["general_plots"]:
        summary_plot_names = crispresso2_info["results"]["general_plots"][
            "summary_plot_names"
        ]
    summary_plot_titles = {}
    if "summary_plot_titles" in crispresso2_info["results"]["general_plots"]:
        summary_plot_titles = crispresso2_info["results"]["general_plots"][
            "summary_plot_titles"
        ]
    summary_plot_labels = {}
    if "summary_plot_labels" in crispresso2_info["results"]["general_plots"]:
        summary_plot_labels = crispresso2_info["results"]["general_plots"][
            "summary_plot_labels"
        ]
    summary_plot_datas = {}
    if "summary_plot_datas" in crispresso2_info["results"]["general_plots"]:
        summary_plot_datas = crispresso2_info["results"]["general_plots"][
            "summary_plot_datas"
        ]

    for plot_name in summary_plot_names:
        if plot_name not in summary_plot_datas:
            summary_plot_datas[plot_name] = {}

    summary_stats = crispresso2_info["results"]["summary_stats"]
    targets = list(summary_stats["target_name"])
    chroms = list(summary_stats["target_chr"])
    positions = list(summary_stats["target_pos"])
    if "region_anno" in summary_stats:
        annotations = list(summary_stats.get("region_anno"))
    else:
        annotations = ["" for _ in targets]
    summary_plot_names = ["targets"] + summary_plot_names
    summary_plot_titles["targets"] = "Targets"
    summary_plot_labels["targets"] = {}
    summary_plot_datas["targets"] = [
        (target, chrom, pos, annotation)
        for target, chrom, pos, annotation in zip(
            targets, chroms, positions, annotations
        )
    ]

    summary_plot_names = ["samples"] + summary_plot_names
    summary_plot_titles["samples"] = "samples"
    summary_plot_labels["samples"] = {}
    summary_plot_datas["samples"] = {"names": {}, "groups": {}}
    for sample in sea_names:
        summary_plot_datas["samples"]["names"][sample] = display_names[sample]
        if groups:
            summary_plot_datas["samples"]["groups"][sample] = groups[sample]

    summary_plot_htmls = {}
    crispresso_data_path = os.path.relpath(
        sea_folder, os.path.dirname(crispressoSea_report_file)
    )
    if crispresso_data_path == ".":
        crispresso_data_path = ""
    else:
        crispresso_data_path += "/"

    sub_html_files = {}
    run_names = []
    for name in sea_names:
        display_name = display_names[name]
        sub_folder = os.path.join("CRISPResso_output", "CRISPRessoPooled_on_" + name)
        crispresso_folder = os.path.join(sea_folder, sub_folder)
        run_data = CRISPRessoShared.load_crispresso_info(
            crispresso_info_file_name=os.path.join(
                crispresso_folder, "CRISPResso2Pooled_info.json"
            )
        )
        if "running_info" not in run_data:
            raise Exception(
                f"CRISPResso run {sub_folder} has no report. Cannot add to Sea report."
            )

        this_sub_html_file = sub_folder + ".html"
        if (
            run_data["running_info"]["args"]
            and run_data["running_info"]["args"].place_report_in_output_folder
        ):
            this_sub_html_file = os.path.join(
                sub_folder, run_data["running_info"]["report_filename"]
            )
        sub_html_files[display_name] = this_sub_html_file

        run_names.append(display_name)

    output_title = "CRISPResso Sea Output"
    if ((crispresso2_info["running_info"]["args"]) and
            ('name' in crispresso2_info["running_info"]["args"]) and
            (crispresso2_info["running_info"]["args"].name != "")):
        output_title += f"<br/>{crispresso2_info['running_info']['args'].name}"

    make_multi_report(
        run_names,
        failed_runs,
        failed_runs_desc,
        sub_html_files,
        crispressoSea_report_file,
        sea_folder,
        _root,
        output_title,
        "sea",
        logger,
        summary_plots={
            "names": summary_plot_names,
            "titles": summary_plot_titles,
            "labels": summary_plot_labels,
            "datas": summary_plot_datas,
            "htmls": summary_plot_htmls,
        },
    )


def reverse_complement(seq):
    """
    Get the reverse complement of a sequence

    params:
    - seq: the sequence to reverse complement

    returns:
    - the reverse complement of the sequence
    """
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n",
    }
    return "".join(complement[base] for base in reversed(seq))


def merge_locations(
    crispresso_pooled_genome_folder,
    genome_file,
    allow_unplaced_chrs=False,
    top_percent_cutoff = 0.2,
    min_amplicon_len=50,
    debug=False,
):
    """
    Merge locations from a CRISPRessoPooled run to get a list of non-overlapping regions that could represent the original amplicons

    params:
    - crispresso_pooled_genome_folder: the folder containing the CRISPRessoPooled output
    - genome_file: path to the genome file
    - allow_unplaced_chrs: whether to allow regions on unplaced chromosomes (chrUn, random, etc). If true, these regions will be included.
    - top_percent_cutoff: the top percent of aligned regions (by region read depth) to consider in finding non-overlapping regions. This is a float between 0 and 1. For example, if set to 0.2, the top 20% of regions (by read depth) will be considered.
    - min_amplicon_len: the minimum length of an amplicon to consider in finding non-overlapping regions. Amplicons shorter than this will be ignored.

    returns:
    - good_region_file: file containing df of merged regions with columns for: chr, start, end, read_count and seq)

    """
    # if output is already completed, just read that in and return it
    good_region_file = os.path.join(crispresso_pooled_genome_folder, "good_nonoverlapping_regions.txt")
    if os.path.isfile(good_region_file):
        info('Using regions from ' + good_region_file, {'percent_complete': 20})
        return good_region_file
    logger.debug("Couldn't find completed regions file at " + good_region_file + ". Generating now.", {'percent_complete': 10})

    #otherwise process the read counts for each region
    genome_read_counts_file = os.path.join(
        crispresso_pooled_genome_folder, "REPORT_READS_ALIGNED_TO_GENOME_ALL_DEPTHS.txt"
    )
    if not os.path.isfile(genome_read_counts_file):
        raise Exception(
            "Genome read counts file not found at " + genome_read_counts_file
        )

    logger.debug('Reading regions from ' + genome_read_counts_file)

    region_info = {}
    good_nonoverlapping_regions = []
    region_df = pd.read_csv(genome_read_counts_file, sep="\t")
    region_df_sub = region_df.sort_values(by="number of reads", ascending=False)

    #filter first by amplicon length
    region_df_sub['length'] = region_df_sub['end'] - region_df_sub['start']
    region_df_sub = region_df_sub[region_df_sub['length'] >= min_amplicon_len]

    top_percent_count = int(len(region_df_sub) * top_percent_cutoff)
    region_df_sub = region_df_sub.head(top_percent_count)
    info(
        "Considering top " + str(top_percent_cutoff * 100) + "% of regions (N=" + str(top_percent_count) + "/" + str(region_df.shape[0]) + ") with at least " 
        + str(region_df_sub["number of reads"].min()) + " reads (mean reads is " + str(region_df_sub["number of reads"].mean()) + ")"
    )

    seen_region_seqs = {}

    for index, region_row in region_df_sub.iterrows():
        region_chr = region_row.chr_id
        region_start = region_row.start
        region_end = region_row.end
        region_count = region_row["number of reads"]
        region_name = "_".join([region_chr, str(region_start), str(region_end)])
        region_seq_output = subprocess.check_output(
            "%s faidx %s %s:%d-%d"
            % (
                "samtools",
                genome_file + ".fa",
                region_chr,
                region_start,
                region_end - 1,
            ),
            shell=True,
        ).decode(sys.stdout.encoding)
        region_seq = "".join(region_seq_output.split("\n")[1:])

        region_info[region_name] = {
            "chr": region_chr,
            "start": region_start,
            "end": region_end,
            "region_count": region_count,
            "seq": region_seq,
        }

        if not allow_unplaced_chrs:
            if "chrUn" in region_chr or "random" in region_chr:
                continue
        region_len = region_end - region_start
        if region_len < min_amplicon_len:
            continue

        if (
            region_seq.upper() in seen_region_seqs
            or reverse_complement(region_seq.upper()) in seen_region_seqs
        ):
            continue

        overlapped_any_region = False
        for other_region in good_nonoverlapping_regions:
            other_region_chr = region_info[other_region]["chr"]
            other_region_start = region_info[other_region]["start"]
            other_region_end = region_info[other_region]["end"]
            overlaps_this_region = False
            if region_chr == other_region_chr:
                if other_region_start <= region_start <= other_region_end:
                    overlaps_this_region = True
                elif other_region_start <= region_end <= other_region_end:
                    overlaps_this_region = True
                elif (
                    other_region_start >= region_start
                    and other_region_end <= region_end
                ):
                    overlaps_this_region = True
                elif (
                    other_region_start <= region_start
                    and other_region_end >= region_end
                ):
                    overlaps_this_region = True
            if overlaps_this_region:
                overlapped_any_region = True
                other_region_count = region_info[other_region]["region_count"]
                if region_count > other_region_count:
                    good_nonoverlapping_regions.remove(other_region)
                    good_nonoverlapping_regions.append(region_name)
        if not overlapped_any_region:
            good_nonoverlapping_regions.append(region_name)
            seen_region_seqs[region_seq.upper()] = True

    total_region_count = region_df.shape[0]
    kept_region_count = len(good_nonoverlapping_regions)
    info("Kept " + str(kept_region_count) + "/" + str(total_region_count) + " regions", {'percent_complete': 20})
    if debug:
        for region in good_nonoverlapping_regions:
            logger.debug(region)
    good_nonoverlapping_region_rows = []
    for region in good_nonoverlapping_regions:
        this_region_info = region_info[region]
        good_nonoverlapping_region_rows.append(
            [
                this_region_info["chr"],
                this_region_info["start"],
                this_region_info["end"],
                this_region_info["region_count"],
                this_region_info["seq"],
            ]
        )

    good_nonoverlapping_region_df = pd.DataFrame(
        good_nonoverlapping_region_rows,
        columns=["chr", "start", "end", "read_count", "seq"],
    )

    info('Wrote good regions to ' + good_region_file, {'percent_complete': 20})
    good_nonoverlapping_region_df.to_csv(good_region_file, sep="\t", index=False)

    return good_region_file

def run_initial_demux(
    experiment_name,
    fastq_r1,
    fastq_r2,
    genome_file,
    output_folder="",
    n_processors=8,
    crispresso_quantification_window_center=-3,
    crispresso_quantification_window_size=1,
    crispresso_base_editor_output=False,
    crispresso_default_min_aln_score=20,
    crispresso_plot_window_size=20,
    allow_unplaced_chrs=False,
    top_percent_cutoff=0.2,
    min_amplicon_len=50,
    suppress_output=True,
    fail_on_pooled_fail=False,
):
    """
    Run CRISPResso on input reads to determine which regions are frequently aligned to

    params:
    - experiment_name: a name for the experiment
    - fastq_r1: path to the R1 fastq file
    - fastq_r2: path to the R2 fastq file
    - genome_file: path to the genome file
    - output_folder: path to the output folder where results should be written
    - n_processors: number of processors to use
    - crispresso_quantification_window_center: the window to use for CRISPResso2 quantification. This will be used for CRISPResso2 subruns and corresponds to the parameter --quantification_window_center. If None, the default value of -3 will be used.
    - crispresso_quantification_window_size: the window size to use for CRISPResso2 quantification. This will be used for CRISPResso2 subruns and corresponds to the parameter --quantification_window_size. If None, the default value of 1 will be used.
    - crispresso_base_editor_output: whether to output base editor plots. If True, base editor plots will be output. If False, base editor plots will not be output.
    - crispresso_default_min_aln_score: the minimum alignment score to use for CRISPResso2. This will be used for CRISPResso2 subruns and corresponds to the parameter --default_min_aln_score. If None, the default value of 20 will be used.
    - crispresso_plot_window_size: the window size to use for CRISPResso2 plots. This will be used for CRISPResso2 subruns and corresponds to the parameter --plot_window_size. If None, the default value of 20 will be used.
    - allow_unplaced_chrs: whether to allow regions on bad chromosomes (chrUn, random, etc). If true, these regions will be included.
    - top_percent_cutoff: the top percent of aligned regions (by region read depth) to consider in finding non-overlapping regions. This is a float between 0 and 1. For example, if set to 0.2, the top 20% of regions (by read depth) will be considered.
    - min_amplicon_len: the minimum length of an amplicon to consider in finding non-overlapping regions. Amplicons shorter than this will be ignored.
    - fail_on_pooled_fail: if true, fail if any pooled CRISPResso run fails. Otherwise, continue even if sub-CRISPResso commands fail.

    returns:
    - merged_regions_file: file with df of merged regions with columns for: chr, start, end, read_count and seq)
    """

    info ("Running initial demultiplexing to find amplicon locations", {'percent_complete': 5})

    r2_string = ""
    if fastq_r2 is not None:
        r2_string = " -r2 " + fastq_r2
    output_string = ""
    if output_folder != "":
        output_string = " -o " + output_folder
    suppress_output_string = ""
    if suppress_output:
        suppress_output_string = " --verbosity 1 "

    stop_on_fail_string = " --skip_failed "
    if not fail_on_pooled_fail:
        stop_on_fail_string = ""

    CRISPResso_output_folder = (
        output_folder + "CRISPRessoPooled_on_" + experiment_name + "_demux"
    )

    command = (
        "CRISPRessoPooled -x " + genome_file
        + " -r1 " + fastq_r1
        + r2_string
        + output_string
        + " --quantification_window_center " + str(crispresso_quantification_window_center) 
        + " --quantification_window_size " + str(crispresso_quantification_window_size)
        + " --default_min_aln_score " + str(crispresso_default_min_aln_score)
        + " -n " + experiment_name + "_demux -p " + str(n_processors)
        + " --no_rerun --keep_intermediate --suppress_plots"
        + " --plot_window_size " + str(crispresso_plot_window_size)
        + suppress_output_string
        + stop_on_fail_string
    )

    crispresso_run_is_complete = False
    if os.path.exists(CRISPResso_output_folder):
        try:
            crispresso_pooled_info = CRISPRessoShared.load_crispresso_info(
                crispresso_info_file_path=CRISPResso_output_folder
                + "/CRISPResso2Pooled_info.json"
            )
            if (
                "demultiplexing_genome_only_regions"
                in crispresso_pooled_info["running_info"]["finished_steps"]
            ):
                crispresso_run_is_complete = True
                info(
                    "CRISPResso output folder already exists for initial demultiplexing, skipping CRISPResso run", {'percent_complete': 40}
                )
        except:
            pass
    if not crispresso_run_is_complete:
        info("Aligning reads to the genome to find amplicon locations")
        debug("Running command " + str(command))
        subprocess.run(command, shell=True, check=True)
        info('Finished running alignment to find amplicon locations', {'percent_complete': 40})

    merged_regions_file = merge_locations(
        CRISPResso_output_folder, genome_file, allow_unplaced_chrs=allow_unplaced_chrs, top_percent_cutoff=top_percent_cutoff, min_amplicon_len=min_amplicon_len, debug=False
    )

    return merged_regions_file


def parse_sample_file(sample_file):
    """
    Parse a sample file and rename columns as neccessary

    params: 
    - sample_file: path to the sample file

    returns:
    - sample_df: a pandas dataframe with the following columns:
        Name: the name of the sample
        fastq_r1: the path to the R1 fastq file
        fastq_r2: the path to the R2 fastq file
        group: the group of the sample
    """
    # Load sample file
    if sample_file.endswith(".xlsx"):
        sample_df = pd.read_excel(sample_file)
    elif sample_file.endswith(".txt"):
        sample_df = pd.read_csv(sample_file, sep="\t")
    else:
        # Attempt to read as tab-delimited if not xlsx or txt
        try:
            sample_df = pd.read_csv(sample_file, sep="\t")
        except Exception as e:
            raise Exception(f"Could not read sample file {sample_file}: {e}")

    renamed_columns = []
    for col in sample_df.columns:
        if col.lower() == 'fastq_r1':
            renamed_columns.append("fastq_r1")
        elif col.lower() == 'r1':
            renamed_columns.append("fastq_r1")
        elif col.lower() == 'fastq_r2':
            renamed_columns.append("fastq_r2")
        elif col.lower() == 'r2':
            renamed_columns.append("fastq_r2")
        elif col.lower() == 'name':
            renamed_columns.append("Name")
        elif col.lower() == 'group':
            renamed_columns.append("group")
        else:
            renamed_columns.append(col)
        
    sample_df.columns = renamed_columns

    #check for unique 'Name' column
    if "Name" not in sample_df.columns:
        raise Exception(
            'Sample file must have column "Name".\n'
            + "Found columns: "
            + str(list(sample_df.columns))
        )
    if sample_df["Name"].duplicated().any():
        raise Exception(
            'In sample file "' + sample_file + '" each sample must have a unique value for "Name".\n'
            + 'Duplicated "Name" values: '
            + str('"' + '", "'.join(sample_df[sample_df["Name"].duplicated()]["Name"].unique()) + '"')
        )

    # add group column if not present
    if 'group' not in sample_df.columns:
        sample_df['group'] = ''
        
    required_columns = ["Name", "fastq_r1"]
    for col in required_columns:
        if col not in sample_df.columns:
            raise Exception(
                'Sample file must have column "'
                + col
                + '".\n'
                + "Found columns: "
                + str(list(sample_df.columns))
                + "\n"
                + "Expecting columns "
                + str(required_columns)
            )
    return sample_df

def parse_target_info(target_file, sort_based_on_mismatch=False):
    """
    Parse a target file listing targets in tab-separated or xls format. Required columns are:
        Guide: the name of the on-target guide
        Sequence: the sequence of the target site (either on-target or off-target)
        PAM: the PAM at the target site
        #MM: the number of mismatches between the on- and off-target. If this is 0 or NA, this target is considered the on-target
        Locus: the locus of the target site in the genome, in the format chr:pos or chr:+pos for positive strand or chr:-pos for negative strand

    Optional column:
    anno or target: an annotation name for the target site

    params:
    target_file: path to the target file
    sort_based_on_mismatch: if true, targets are sorted based on mismatch count. If false, targets are presented in the order they are in the file

    The On-target for each target is the target with the fewest mismatches (may not necessarily be 0MM for retargeting guides). If there are multiple targets with the same number of mismatches, the first one is used.

    Returns:
    - a pandas dataframe with the following columns for processing:
        target_id: a unique identifier for the target
        target_name: the name of the target
        sort_index: the index of the target in the original file
        target_chr: the chromosome of the target
        target_pos: the position of the target
        target_seq_no_gaps_with_pam: the target sequence without gaps and with PAM
        target_seq_no_gaps: the target sequence without gaps
        target_seq_with_gaps: the target sequence with gaps
        target_pam: the target pam sequence
        ontarget_name: the name of the on-target guide
        ontarget_sequence: the sequence of the on-target guide
    """
    if not os.path.exists(target_file):
        raise Exception("Target file not found at " + target_file)

    if target_file.endswith(".xlsx"):
        target_df = pd.read_excel(target_file)
    else:
        target_df = pd.read_csv(target_file, sep="\t")

    required_columns = ["Guide", "Sequence", "PAM", "#MM", "Locus"]

    renamed_columns = []
    for col in target_df.columns:
        col_is_in_required = False
        for required_col in required_columns:
            if col.lower() == required_col.lower():
                renamed_columns.append(required_col)
                col_is_in_required = True
                break
        if col_is_in_required:
            pass
        elif col.lower() == 'seq':
            renamed_columns.append("Sequence")
        elif col.lower() == 'mismatches' or col.lower() == 'mm':
            renamed_columns.append("#MM")
        elif col.lower() == 'loc' or col.lower() == 'location':
            renamed_columns.append("locus")
        else:
            renamed_columns.append(col)
        
    target_df.columns = renamed_columns

    for col in required_columns:
        if col not in target_df.columns:
            raise Exception(
                'Target file must have column "' + col + '".\n'
                + "Found columns: " + str(list(target_df.columns)) + "\n"
                + "Expecting columns " + str(required_columns)
            )

    target_df["#MM"] = target_df["#MM"].fillna(0).astype(int)
    target_df["sort_index"] = target_df.index

    on_targets = target_df.sort_values(by=["#MM", "sort_index"]).drop_duplicates(
        subset="Guide", keep="first"
    )
    if sort_based_on_mismatch:
        target_df["sort_index"] = list(range(1, target_df.shape[0]))

    for idx, row in on_targets.iterrows():
        if row["#MM"] != 0:
            warn(
                "Warning: On-target target "
                + row["Guide"]
                + " has "
                + str(row["#MM"])
                + " mismatches (expecting 0 mismatches for on-target)"
            )

    ontarget_seqs = {}
    for idx, row in on_targets.iterrows():
        ontarget_seqs[row["Guide"]] = row["Sequence"]

    target_df["ontarget_name"] = "NA"
    target_df["ontarget_sequence"] = "NA"
    target_df["target_id"] = "NA"  # unique id
    target_df["target_name"] = "NA"  # may be supplied by user
    target_df["target_chr"] = "NA"
    target_df["target_pos"] = "NA"
    target_df["target_seq_with_gaps"] = target_df["Sequence"]
    target_df["target_pam"] = target_df["PAM"]
    target_df["target_seq_no_gaps"] = "NA"
    target_df["target_seq_no_gaps_with_pam"] = "NA"
    for idx, row in target_df.iterrows():
        this_ontarget_name = row["Guide"]
        if this_ontarget_name not in ontarget_seqs:
            raise Exception(
                "On-target sequence not found for target "
                + this_ontarget_name
                + ". Found ontargets: "
                + str(list(ontarget_seqs.keys()))
            )
        if str(row["#MM"]) == "0":
            this_target_id = f'{idx}_{this_ontarget_name}_ON'
        else:
            this_target_id = f'{idx}_{this_ontarget_name}_OB{str(int(row["#MM"]))}'

        target_df.loc[idx, "target_id"] = this_target_id

        this_target_name = this_target_id
        if "target" in row and row["target"] is not None and str(row["target"]) != "nan":
            this_target_name = row["target"]
        if "Target" in row and row["Target"] is not None and str(row["Target"]) != "nan":
            this_target_name = row["Target"]
        if "anno" in row and row["anno"] is not None and str(row["anno"]) != "nan":
            this_target_name = row["anno"]
        target_df.loc[idx, "target_name"] = this_target_name

        if this_ontarget_name not in ontarget_seqs:
            raise Exception(
                "Cannot find ontarget name for "
                + str(this_ontarget_name)
                + " in "
                + str(ontarget_seqs.keys())
            )
        this_ontarget_seq = ontarget_seqs[this_ontarget_name]
        target_df.loc[idx, "ontarget_sequence"] = this_ontarget_seq
        target_df.loc[idx, "ontarget_name"] = this_ontarget_name

        this_chr_loc_els = row["Locus"].split(":")
        this_chr = this_chr_loc_els[0]
        if '+' in this_chr_loc_els[1] or '-' in this_chr_loc_els[1]:
            this_pos = int(this_chr_loc_els[1][1:])
        else:
            this_pos = int(this_chr_loc_els[1])
        target_df.loc[idx, "target_chr"] = this_chr
        target_df.loc[idx, "target_pos"] = this_pos
        target_df.loc[idx, "target_seq_no_gaps_with_pam"] = (
            row["Sequence"].replace("-", "") + row["PAM"]
        )
        target_df.loc[idx, "target_seq_no_gaps"] = row["Sequence"].replace("-", "")

    return target_df[
            [
                "target_id",
                "target_name",
                "sort_index",
                "target_chr",
                "target_pos",
                "target_seq_no_gaps_with_pam",
                "target_seq_no_gaps",
                "target_seq_with_gaps",
                "target_pam",
                "ontarget_name",
                "ontarget_sequence",
            ]
        ]


def make_target_region_assignments(
    merged_region_info_file,
    target_file,
    output_folder,
    sort_based_on_mismatch=False,
    allow_target_match_to_other_region_loc=False,
):
    """
    Assign guide targets to regions based on their sequences and positions

    params:
    - merged_region_info_df: dataframe of merged regions with columns for: chr, start, end, read_count and seq
    - target_file: path to the user-provided target file
    - output_folder: path to the output file root
    - sort_based_on_mismatch: if true, targets are sorted based on mismatch count. If false, targets are presented in the order they are in the file
    - allow_target_match_to_other_region_loc: if true, targets can be matched to regions based on sequence even if they are not in the same position. If false, targets can only match to regions containing that target's chr:loc

    returns:
    - path to the CRISPRessoPooledRegions file
    - target_df: the target dataframe
    - region_df: the region dataframe
    """
    # we parse the target_df here because we'll add more columns to it
    target_df = parse_target_info(target_file, sort_based_on_mismatch)
    merged_region_df = pd.read_csv(merged_region_info_file, sep="\t")
    merged_region_df.index = merged_region_df["chr"].astype(str) + "_" + merged_region_df["start"].astype(str) + "_" + merged_region_df["end"].astype(str)

    # first, for each target, match its sequence to genomic region(s)
    # the matching will be performed so that each target is assigned to
    # 1) If there is one region that corresponds to the target position it is assigned to that region.
    # 2) Else if there are multiple regions that correspond to the target position, it is assigned to the region with the highest read count
    # 3) Else (the chr:pos is not given or the target position is in no regions), if the target sequence is found in one region the target is assigned to that region
    # 4) Else (there are multiple regions that contain the target sequence), it is assigned to the region with the highest read count without a previously-assigned target.
    target_matches = {}  # target_id -> final assigned region
    region_matches = defaultdict(
        list
    )  # region -> list of target_ids that were assigned to it

    all_target_matches = defaultdict(
        list
    )  # target_id -> list of regions (all regions with position or sequence match)
    all_region_matches = defaultdict(list)  # region -> list of target_ids

    for target_idx, target_row in target_df.iterrows():
        target_seq_with_pam = target_row["target_seq_no_gaps_with_pam"]
        target_id = target_row["target_id"]
        target_name = target_row["target_name"]
        target_chr = target_row["target_chr"]
        target_pos = target_row["target_pos"]

        this_target_matches = []
        this_target_region_matches = []
        for region_idx, region_row in merged_region_df.iterrows():
            region_name = region_row.name
            region_chr = region_row["chr"]
            region_start = region_row["start"]
            region_end = region_row["end"]
            region_seq = region_row["seq"]

            # first, check if the region contains the target sequence
            is_seq_match = False
            if target_seq_with_pam.upper() in region_seq.upper():
                is_seq_match = True
            elif reverse_complement(target_seq_with_pam).upper() in region_seq.upper():
                is_seq_match = True

            if (
                not allow_target_match_to_other_region_loc
            ):  # don't allow matches based on target match to region sequence
                is_seq_match = False

            if is_seq_match:
                this_target_matches.append(region_name)

            # next, check if the target position is within the region
            is_pos_match = False
            if (
                str(target_chr) == str(region_chr)
                and target_pos >= region_start
                and target_pos <= region_end
            ):
                this_target_region_matches.append(region_name)
                is_pos_match = True
            elif (
                str(target_chr).replace("chr", "") == str(region_chr).replace("chr", "")
                and target_pos >= region_start
                and target_pos <= region_end
            ):
                this_target_region_matches.append(region_name)
                is_pos_match = True

            if is_seq_match or is_pos_match:
                all_target_matches[target_id].append(region_name)
                all_region_matches[region_name].append(target_id)

        # make the final assignment for this target
        if len(this_target_region_matches) == 1:  # if there is one region match, take it
            target_matches[target_id] = this_target_region_matches[0]
            region_matches[this_target_region_matches[0]].append(target_id)
        elif (
            len(this_target_region_matches) > 1
        ):  # if there are multiple region matches, take the one with highest number of reads
            best_region_name = None
            best_region_count = 0
            for region_name in this_target_region_matches:
                region_count = merged_region_df.loc[region_name, "region_count"]
                if region_count > best_region_count:
                    best_region_name = region_name
                    best_region_count = region_count
            target_matches[target_id] = best_region_name
            region_matches[best_region_name].append(target_id)
        elif (
            len(this_target_matches) == 1
        ):  # if there is one sequence match, take it
            target_matches[target_id] = this_target_matches[0]
            region_matches[this_target_matches[0]].append(target_id)
        elif (
            len(this_target_matches) > 1
        ):  # if there are multiple sequence matches, take the one with highest number of reads without a target match
            best_region_name = this_target_matches[0]
            best_region_count = 0
            for region_name in this_target_matches:
                region_count = merged_region_df.loc[region_name, "region_count"]
                if (
                    region_count > best_region_count
                    and len(region_matches[region_name]) == 0
                ):
                    best_region_name = region_name
                    best_region_count = region_count
            target_matches[target_id] = best_region_name
            region_matches[best_region_name].append(target_id)

    # next, finalize matches and write them to file
    out_file = output_folder + "target_region_matches.txt"
    target_match_count = 0  # how many targets had a region match
    target_nomatch_count = 0  # how many targets did not have a region match

    target_df["matched_region_count"] = 0
    target_df["matched_region"] = "NA"
    target_df["matched_region_id"] = "NA"
    with open(out_file, "w") as fout:
        fout.write(
            "target_id\ttarget_name\ttarget_seq\ttarget_seq_with_pam\ttarget_chr\ttarget_pos\tontarget_name\tontarget_seq\tmatched_region_count\tmatched_regions\tregion_chr\tregion_start\tregion_end\tregion_read_count\tregion_seq\n"
        )
        for target_idx, target_row in target_df.iterrows():
            target_id = target_row["target_id"]
            target_seq = target_row["target_seq_no_gaps"]
            target_seq_with_pam = target_row["target_seq_no_gaps_with_pam"]
            target_name = target_row["target_name"]
            target_chr = target_row["target_chr"]
            target_pos = target_row["target_pos"]
            target_ontarget_name = target_row["ontarget_name"]
            target_ontarget_seq = target_row["ontarget_sequence"]
            if target_id in target_matches:
                target_match_count += 1
                matched_region_count = len(all_target_matches[target_id])
                matched_regions = ", ".join(all_target_matches[target_id])
                region_name = target_matches[target_id]
                region_chr = merged_region_df.loc[region_name, "chr"]
                region_start = merged_region_df.loc[region_name, "start"]
                region_end = merged_region_df.loc[region_name, "end"]
                region_count = merged_region_df.loc[region_name, "read_count"]
                region_seq = merged_region_df.loc[region_name, "seq"]
                fout.write(
                    "\t".join(
                        [
                            str(x)
                            for x in [
                                target_id,
                                target_name,
                                target_seq,
                                target_seq_with_pam,
                                target_chr,
                                target_pos,
                                target_ontarget_name,
                                target_ontarget_seq,
                                matched_region_count,
                                matched_regions,
                                region_chr,
                                region_start,
                                region_end,
                                region_count,
                                region_seq,
                            ]
                        ]
                    )
                    + "\n"
                )
                target_df.loc[target_idx, "matched_region_count"] = matched_region_count
                target_df.loc[target_idx, "matched_region"] = region_name
            # if target wasn't selected as the final match for a region, it will be in all_target_matches
            elif target_id in all_target_matches:
                target_nomatch_count += 1
                matched_region_count = len(all_target_matches[target_id])
                matched_regions = ", ".join(all_target_matches[target_id])
                fout.write(
                    "\t".join(
                        [
                            str(x)
                            for x in [
                                target_id,
                                target_name,
                                target_seq,
                                target_seq_with_pam,
                                target_chr,
                                target_pos,
                                target_ontarget_name,
                                target_ontarget_seq,
                                0,
                                matched_regions,
                                "NA",
                                "NA",
                                "NA",
                                "NA",
                                "NA",
                            ]
                        ]
                    )
                    + "\n"
                )
            else:
                target_nomatch_count += 1
                fout.write(
                    "\t".join(
                        [
                            str(x)
                            for x in [
                                target_id,
                                target_name,
                                target_seq,
                                target_seq_with_pam,
                                target_chr,
                                target_pos,
                                target_ontarget_name,
                                target_ontarget_seq,
                                0,
                                "NA",
                                "NA",
                                "NA",
                                "NA",
                                "NA",
                                "NA",
                            ]
                        ]
                    )
                    + "\n"
                )

    info("Matched " + str(target_match_count) + "/" + str(len(target_df)) + " targets to regions from read alignments. Wrote matches to " + out_file, {'percent_complete': 25})

    # next write information about all regions plus the file for crispresso input
    crispresso_output_file = output_folder + "CRISPRessoPooledRegions.txt"
    all_region_output_file = output_folder + "all_regions.txt"
    region_match_count = 0
    region_nomatch_count = 0
    printed_crispresso_count = 0
    with open(all_region_output_file, "w") as rout, open(crispresso_output_file, "w") as cout:
        rout.write(
            "region_id\tchr\tstart\tend\tread_count\tseq\ttarget_match_count\ttarget_matches\ttarget_id\ttarget_name\ttarget_seq\ttarget_chr\ttarget_pos\n"
        )
        for region_idx, region_row in merged_region_df.iterrows():
            region_chr = region_row["chr"]
            region_start = region_row["start"]
            region_end = region_row["end"]
            region_seq = region_row["seq"]
            region_count = region_row["read_count"]
            region_name = region_row.name

            target_id = "NA"
            target_name = "NA"
            target_seq = "NA"
            target_chr = "NA"
            target_pos = "NA"
            this_target_match_count = len(all_region_matches[region_name])
            if this_target_match_count > 0:
                this_target_matches = ", ".join(all_region_matches[region_name])
            else:
                this_target_matches = "NA"
            if region_name in all_region_matches:
                if region_name in region_matches:
                    region_match_count += 1
                    target_id = region_matches[region_name][0]  # the assigned targets that matched to this region
                    target_seq = target_df.loc[target_df["target_id"] == target_id, "target_seq_no_gaps"].values[0]  # no pam (for CRISPResso)
                    target_name = target_df.loc[target_df["target_id"] == target_id, "target_name"].values[0]
                    target_chr = target_df.loc[target_df["target_id"] == target_id, "target_chr"].values[0]
                    target_pos = target_df.loc[target_df["target_id"] == target_id, "target_pos"].values[0]
                    region_name = region_name + "_" + target_name
                    target_df.loc[target_df["target_id"] == target_id, "matched_region_id"] = region_name
                    cout.write("\t".join([region_name, region_seq, target_seq]) + "\n")
                    printed_crispresso_count += 1
                else:
                    region_nomatch_count += 1
            else:
                region_nomatch_count += 1
            rout.write(
                "\t".join(
                    str(x)
                    for x in [ region_name, region_chr, region_start, region_end, region_count, region_seq, 
                              this_target_match_count, this_target_matches, target_id, target_name, 
                              target_seq, target_chr, target_pos, ]
                )
                + "\n"
            )

    info("Matched " + str(region_match_count) + "/" + str(region_match_count + region_nomatch_count) + " frequently-aligned locations to targets. " + \
        "Wrote region info to " + all_region_output_file, {'percent_complete': 30})

    if printed_crispresso_count == 0:
        raise Exception(
            "No regions matched to targets. Check the target file and the regions file for consistency.\n" + \
            "\tTarget file: " + target_file + "\n\tRegions file: " + merged_region_info_file
        )
    region_df = pd.read_csv(all_region_output_file, sep="\t")
    return crispresso_output_file, target_df, region_df

def run_crispresso_with_assigned_regions(
    experiment_name,
    fastq_r1,
    fastq_r2,
    genome_file,
    crispresso_region_file,
    output_folder,
    n_processors=8,
    crispresso_quantification_window_center=-3,
    crispresso_quantification_window_size=1,
    crispresso_base_editor_output=False,
    crispresso_default_min_aln_score=60,
    crispresso_plot_window_size=20,
    suppress_commandline_output=True,
    fail_on_pooled_fail=False,
):
    """
    Run CRISPResso on a specific sample with the assigned regions

    params:
    - experiment_name: a name for the experiment
    - fastq_r1: path to the R1 fastq file
    - fastq_r2: path to the R2 fastq file
    - genome_file: path to the genome file
    - crispresso_region_file: path to the file with region assignments
    - n_processors: number of processors to use
    - crispresso_quantification_window_center: the center of the quantification window
    - crispresso_quantification_window_size: the size of the quantification window
    - crispresso_base_editor_output: whether to output base editor results
    - crispresso_default_min_aln_score: the minimum alignment score for CRISPResso
    - crispresso_plot_window_size: the size of the plot window
    - suppress_commandline_output: whether to suppress commandline output from CRISPRessoPooled
    - fail_on_pooled_fail: if true, fail if any pooled CRISPResso run fails. Otherwise, continue even if sub-CRISPResso commands fail.

    returns:
    - the output folder for this CRISPResso run
    """
    # Run CRISPResso again on the subset, this time including target sequences
    r2_string = ""
    if fastq_r2 is not None:
        r2_string = " -r2 " + fastq_r2
    output_string = ""
    if output_folder != "":
        output_string = " -o " + output_folder

    suppress_output_string = ""
    if suppress_commandline_output:
        suppress_output_string = " --verbosity 1"

    stop_on_fail_string = ""
    if not fail_on_pooled_fail:
        stop_on_fail_string = " --skip_failed"

    base_editor_string = ""
    if crispresso_base_editor_output:
        base_editor_string = " --base_editor_output "

    CRISPResso_output_folder = output_folder + "CRISPRessoPooled_on_" + experiment_name

    command = (
        "CRISPRessoPooled -f "
        + crispresso_region_file
        + " -x " + genome_file
        + " -n " + experiment_name
        + " -r1 " + fastq_r1
        + r2_string
        + output_string
        + " --min_reads_to_use_region 1 "
        + base_editor_string 
        + " --quantification_window_center " + str(crispresso_quantification_window_center) 
        + " --quantification_window_size " + str(crispresso_quantification_window_size)
        + " --default_min_aln_score " + str(crispresso_default_min_aln_score)
        + " -p " + str(n_processors)
        + " --no_rerun  --exclude_bp_from_left 0 --exclude_bp_from_right 0 "
        + " --plot_window_size " + str(crispresso_plot_window_size)
        + suppress_output_string
        + stop_on_fail_string
    )

    crispresso_run_is_complete = False
    if os.path.exists(CRISPResso_output_folder):
        try:
            crispresso_pooled_info = CRISPRessoShared.load_crispresso_info(
                crispresso_info_file_path=CRISPResso_output_folder + "/CRISPResso2Pooled_info.json"
            )
            if "end_time" in crispresso_pooled_info["running_info"]:
                crispresso_run_is_complete = True
                info("CRISPResso output already exists for sample " + experiment_name + ", skipping CRISPResso run")
        except:
            raise Exception("failed!")
            pass

    if not crispresso_run_is_complete:
        info("Aligning reads to the genome to find amplicon locations")
        info("Running command " + str(command))
        try:
            subprocess.run(command, shell=True, check=True)
        except:
            info("Failed to run CRISPRessoPooled for sample " + experiment_name)
            return None

    return CRISPResso_output_folder


def analyze_run(output_name, crispresso_folder, target_df, region_df):
    """
    Analyzes a CRISPResso Pooled run to get the base editing and indel rates for each site

    params:
    - output_name: the name of the output file
    - crispresso_folder: the folder containing the CRISPResso output
    - target_df: the target dataframe with columns target_id, target_name, matched_region, matched_region_id, etc.
    - region_df: the region dataframe with columns region_id, target_id, etc.

    returns:
    - summary_output_file: a file with the summary of the run
    - target_summary_df: a pandas data frame with results for each target
    """
    target_summary_file = output_name + ".complete_target_summary.txt"

    # if no output folder is provided, add row with NAs for each target
    if crispresso_folder is None or not os.path.exists(crispresso_folder):
        target_data = []
        with open(target_summary_file, "w") as fout:
            fout.write(
                "target_id\ttarget_label\tpooled_result_name\thighest_a_g_pct\thighest_c_t_pct\thighest_indel_pct\ttot_reads\n"
            )
        for idx, row in target_df.iterrows():
            target_id = row["target_id"]
            target_name = row["target_name"]
            fout.write(target_id + "\t" + target_name + "\tNA\tNA\tNA\tNA\n")
            target_data.append([target_id, target_name, np.nan, np.nan, np.nan, np.nan])

        target_data_df = pd.DataFrame(
            target_data,
            columns=[
                "target_id",
                "target_label",
                "pooled_result_name",
                "highest_a_g_pct",
                "highest_c_t_pct",
                "highest_indel_pct",
                "tot_reads",
            ],
        )
        return target_summary_file, target_data_df

    pooled_info = CRISPRessoShared.load_crispresso_info(
        crispresso_info_file_path=crispresso_folder + "/CRISPResso2Pooled_info.json"
    )
    pooled_results = pooled_info["results"]["good_region_folders"] # dict of result_name to CRISPResso folder

    if len(pooled_results) == 0:
        raise Exception(
            f"No results found in CRISPResso output folder {crispresso_folder}. Check the output folder for errors."
        )

    output_dir = output_name + "_plots"
    os.makedirs(output_dir, exist_ok=True)

    #associate pooled results with targets
    target_pooled_results_lookup = {}
    for pooled_result_name in pooled_results.keys():
        folder_name = pooled_results[pooled_result_name]
        run_info = CRISPRessoShared.load_crispresso_info(
            crispresso_info_file_path=crispresso_folder
            + "/"
            + folder_name
            + "/CRISPResso2_info.json"
        )
        region_info = region_df.loc[region_df["region_id"] == pooled_result_name]
        target_id = region_info["target_id"].values[0]

        target_pooled_results_lookup[target_id] = pooled_result_name

    #iterate through targets in order and aggregate results
    target_data = []
    with open(target_summary_file, "w") as fout:
        fout.write(
            "target_id\ttarget_label\tpooled_result_name\thighest_a_g_pct\thighest_c_t_pct\thighest_indel_pct\ttot_reads\n"
        )
        for idx, row in target_df.iterrows():
            target_id = row["target_id"]
            target_name = row["target_name"]
            if target_id not in target_pooled_results_lookup:
                fout.write(target_id + "\t" + target_name + "\tNA\tNA\tNA\tNA\n")
                target_data.append(
                    [target_id, target_name, np.nan, np.nan, np.nan, np.nan]
                )
            else:
                pooled_result_name = target_pooled_results_lookup[target_id]
                folder_name = pooled_results[pooled_result_name]

                run_info = CRISPRessoShared.load_crispresso_info(
                    crispresso_info_file_path=crispresso_folder
                    + "/"
                    + folder_name
                    + "/CRISPResso2_info.json"
                )
                region_info = region_df.loc[region_df["region_id"] == pooled_result_name]

                target_info = target_df.loc[target_df["target_id"] == target_id]
                target_name = target_info["target_name"].values[0]

                include_idxs = run_info["results"]["refs"]["Reference"]["include_idxs"]
                nuc_freq_table_file = (crispresso_folder + "/" + folder_name + "/Nucleotide_frequency_table.txt")
                mod_freq_table_file = (crispresso_folder + "/" + folder_name + "/Modification_count_vectors.txt")
                nuc_freq_table = pd.read_csv(nuc_freq_table_file, sep="\t", index_col=0)
                mod_freq_table = pd.read_csv(mod_freq_table_file, sep="\t", index_col=0)
                highest_a_g_pct = 0
                has_valid_a_g = False
                highest_c_t_pct = 0
                has_valid_c_t = False
                highest_indel_pct = 0
                tot_count = 0
                for this_idx in include_idxs:
                    nuc_vals = nuc_freq_table.iloc[:, this_idx]
                    mod_vals = mod_freq_table.iloc[:, this_idx]

                    this_nuc = nuc_vals.name # e.g. C.15 -- the 15 comes from pandas making column names unique
                    a_count = nuc_vals["A"]
                    c_count = nuc_vals["C"]
                    g_count = nuc_vals["G"]
                    t_count = nuc_vals["T"]
                    del_count = mod_vals["Deletions"]
                    ins_count = mod_vals["Insertions"]
                    tot_count = mod_vals["Total"]

                    if this_nuc[0].upper() == "A":
                        has_valid_a_g = True
                        a_g_pct = g_count / tot_count
                        if a_g_pct > highest_a_g_pct:
                            highest_a_g_pct = a_g_pct
                    elif this_nuc[0].upper() == "C":
                        has_valid_c_t = True
                        c_t_pct = t_count / tot_count
                        if c_t_pct > highest_c_t_pct:
                            highest_c_t_pct = c_t_pct
                    elif this_nuc[0].upper() == "G":
                        has_valid_c_t = True
                        c_t_pct = a_count / tot_count
                        if c_t_pct > highest_c_t_pct:
                            highest_c_t_pct = c_t_pct
                    elif this_nuc[0].upper() == "T":
                        has_valid_a_g = True
                        a_g_pct = c_count / tot_count
                        if a_g_pct > highest_a_g_pct:
                            highest_a_g_pct = a_g_pct
                    if (del_count + ins_count) / tot_count > highest_indel_pct:
                        highest_indel_pct = (del_count + ins_count) / tot_count

                highest_a_g_pct *= 100
                highest_c_t_pct *= 100
                highest_indel_pct *= 100

                if not has_valid_a_g:
                    highest_a_g_pct = np.nan
                if not has_valid_c_t:
                    highest_c_t_pct = np.nan

                data_to_write = [
                    target_id,
                    target_name,
                    pooled_result_name,
                    highest_a_g_pct,
                    highest_c_t_pct,
                    highest_indel_pct,
                    tot_count,
                ]
                fout.write("\t".join(str(x) for x in data_to_write) + "\n")
                target_data.append(data_to_write)

            #copy over CRISPResso plot 2b for this target
            plot2b_file = (crispresso_folder + "/" + folder_name + "/" + run_info["results"]["refs"]["Reference"]["plot_2b_roots"][0] + ".pdf")
            if not os.path.isfile(plot2b_file):
                warn("WARNING: missing plot2b file for", pooled_result_name)
                continue
            destination_name = target_name + "." + pooled_result_name + ".plot2b.pdf"
            os.popen("cp " + plot2b_file + " " + output_dir + "/" + destination_name)

    plot_heatmap(output_name)
    target_data_df = pd.DataFrame(
        target_data,
        columns=[
            "target_id",
            "target_label",
            "pooled_result_name",
            "highest_a_g_pct",
            "highest_c_t_pct",
            "highest_indel_pct",
            "tot_reads",
        ],
    )
    return target_summary_file, target_data_df


def view_complete_target_summary(output_name):
    """
    View the complete_target_summary.txt file for debugging

    params:
    - output_name: the name of the output folder

    """
    d = pd.read_csv(output_name + ".complete_target_summary.txt", sep="\t")
    d.dropna(inplace=True)
    d.index = d["target_id"]
    print(d)


def plot_heatmap(output_name):
    """
    Plot a heatmap of the highest_a_g_pct, highest_c_t_pct, and highest_indel_pct for each target
    Data is read from output_name + ".complete_target_summary.txt"
    Plot is written to output_name + ".heatmap.pdf"

    params:
    - output_name: the name of the output folder
    """
    d = pd.read_csv(output_name + ".complete_target_summary.txt", sep="\t")
    d.dropna(inplace=True)
    d.set_index("target_id", inplace=True)
    fig = plt.figure(figsize=(12, 12), dpi=100, facecolor="w", edgecolor="k")
    sns.heatmap(
        d[["highest_a_g_pct", "highest_c_t_pct", "highest_indel_pct"]].astype(float),
        annot=True,
        fmt=".4f",
        cmap="Reds",
    )
    fig.savefig(output_name + ".heatmap.pdf")
    plt.close(fig)


def plot_heatmap_sub(output_name, target_name):
    """
    Plot a heatmap of the highest_a_g_pct, highest_c_t_pct, and highest_indel_pct for a specific target
    Data is read from output_name + ".complete_target_summary.txt"
    Plot is written to output_name + ".heatmap.pdf"

    params:
    - output_name: the name of the output folder
    - target_name: the name of the target to plot
    """
    d = pd.read_csv(output_name + ".complete_target_summary.txt", sep="\t")
    d.dropna(inplace=True)
    d.set_index("target_id", inplace=True)
    dsub = d[d["ontarget_name"] == target_name]
    fig = plt.figure(figsize=(6, 6), dpi=100, facecolor="w", edgecolor="k")
    sns.heatmap(
        dsub[["highest_a_g_pct", "highest_c_t_pct", "highest_indel_pct"]].astype(float),
        annot=True,
        fmt=".4f",
        cmap="Reds",
    )
    plt.close(fig)


def create_target_df_for_plotting(target_df):
    """
    Create a dataframe for plotting target in a heatmap

    params:
    - target_df: the target dataframe
        includes columns: target_id, target_name, target_seq_no_gaps, target_pam, ontarget_name, ontarget_sequence

    returns:
    - a dataframe with target sequences for plotting
        adds column: padded_seqs_for_plot (list of characters for each box in the plot heatmap)
    """
    aln_matrix = CRISPResso2Align.make_matrix()

    seqs_for_plot = []
    aln_ontarget_seqs = []
    aln_offtarget_seqs = []

    last_target_name = None # target name - not the on/offtarget name
    last_target_seq = None
    for idx, row in target_df.iterrows():
        if row["ontarget_name"] != last_target_name: # this is the first time we see this target sequence - plot the full sequence
            last_target_name = row["ontarget_name"]
            last_target_seq = row["target_seq_no_gaps"] # may change this to row['ontarget_sequence'] if everything should be plotted with regard to the on-target sequence. As is, everything is plotted with regard to the vertically first time we see the target sequence
            this_seq = row["target_seq_no_gaps"]
            seqs_for_plot.append((list(this_seq), row["target_pam"]))
            aln_ontarget_seqs.append(this_seq)
            aln_offtarget_seqs.append(this_seq)
        else: #otherwise, only plot only '.' or differences compared to the first time we see this target sequence
            ref_incentive = np.zeros(len(last_target_seq) + 1, dtype=int)
            off_target_aln, on_target_aln, aln_score = CRISPResso2Align.global_align(
                row["target_seq_no_gaps"].upper(),
                last_target_seq.upper(),
                matrix=aln_matrix,
                gap_incentive=ref_incentive,
                gap_open=-20,
                gap_extend=-20,
            )

            this_chars_for_plot = []
            first_base_is_del_in_ontarget = (
                False  # toggle if starts with run of deletions in ontarget
            )
            # bases are collapsed left meaning that an insertion will be grouped with the base to its left (except for the case when a deletion starts, and it is grouped to the right)
            # example:
            #  on_target_aln  = "---A-A-CCTG--"
            #  off_target_aln = "ACTAGATCCTGCC"
            # produces:
            # ACTAG,AT,.,.,.,GCC
            for idx, (off_target_char, on_target_char) in enumerate(
                zip(off_target_aln, on_target_aln)
            ):
                if on_target_char == off_target_char:
                    if (
                        first_base_is_del_in_ontarget
                    ):  # if we finished a run of starting deletions, add this base to the previous char
                        this_chars_for_plot[-1] += off_target_char
                        first_base_is_del_in_ontarget = False
                    else:
                        this_chars_for_plot.append(".")
                elif (
                    on_target_char == "-"
                ):  # if gap in ontarget (insertion in this target/offtarget), add this base to the previous base
                    if idx == 0:
                        first_base_is_del_in_ontarget = True
                        this_chars_for_plot.append(off_target_char)
                    else:
                        if this_chars_for_plot[-1] == ".":
                            this_chars_for_plot[-1] = on_target_aln[idx - 1]
                        this_chars_for_plot[-1] += off_target_char
                else:  # mismatch
                    if (
                        first_base_is_del_in_ontarget
                    ):  # if we finished a run of starting deletions, add this base to the previous char
                        this_chars_for_plot[-1] += off_target_char
                        first_base_is_del_in_ontarget = False
                    else:
                        this_chars_for_plot.append(off_target_char)

            post_len = len("".join(this_chars_for_plot))
            if post_len != len(off_target_aln):
                info(
                    "offtarget lengths do not match ("
                    + str(post_len)
                    + " vs "
                    + str(len(off_target_aln))
                    + ")"
                )
            if len(this_chars_for_plot) != len(on_target_aln.replace("-", "")):
                info(
                    "ontarget lengths do not match ("
                    + str(len(this_chars_for_plot))
                    + " vs "
                    + str(len(on_target_aln.replace("-", "")) + ")")
                )

            seqs_for_plot.append((this_chars_for_plot, row["target_pam"]))
            aln_ontarget_seqs.append(on_target_aln)
            aln_offtarget_seqs.append(off_target_aln)

    max_target_seq_len = 0
    max_pam_len = 0
    for target_seq_chars, pam_seq in seqs_for_plot:
        if len(target_seq_chars) > max_target_seq_len:
            max_target_seq_len = len(target_seq_chars)
        if len(pam_seq) > max_pam_len:
            max_pam_len = len(pam_seq)

    # now, pad the alignments with spaces in case they are different lengths
    padded_seqs_for_plot = []
    for target_seq_chars, pam_seq in seqs_for_plot:
        padded_target_chars = [" "] * (max_target_seq_len - len(target_seq_chars)) + target_seq_chars
        padded_pam_seq = pam_seq + " " * (max_pam_len - len(pam_seq))
        padded_seqs_for_plot.append(padded_target_chars + list(padded_pam_seq))

    target_df["padded_seqs_for_plot"] = padded_seqs_for_plot
    # these are for debugging - can be displayed as a table
    # target_df['aln_ontarget_seqs'] = aln_ontarget_seqs
    # target_df['aln_offtarget_seqs'] = aln_offtarget_seqs

    df_targets = pd.DataFrame(
        target_df["padded_seqs_for_plot"]
        .apply(lambda seq: [char for char in seq])
        .tolist()
    )
    df_targets.index = target_df["target_id"]

    return df_targets


def create_plots(
    data_df,
    sample_df,
    target_plot_df,
    output_folder=None,
    file_prefix="CRISPRessoSea",
    crispresso2_info=None,
    heatmap_fig_height=24,
    heatmap_fig_width=24,
    heatmap_seq_plot_ratio=1,
    heatmap_title_fontsize=30,
    heatmap_y_tick_fontsize=16,
    heatmap_x_tick_fontsize=16,
    heatmap_nucleotide_fontsize=14,
    heatmap_legend_fontsize=20,
    heatmap_legend_ncol=None,
    heatmap_max_value=None,
    heatmap_min_value=None,
    plot_group_order=None,
    dot_plot_ylims=[None, None],
    sig_method_parameters=None,
):
    """For each group, plot the highest_a_g_pct, highest_c_t_pct, highest_indel_pct, and tot_reads for each target

    Args:
        data_df (pd.DataFrame): a dataframe with data for plotting
            includes columns: target_id, target_name, target_seq_no_gaps_with_pam, target_chr, target_pos, region_anno
        sample_df (pd.DataFrame): a dataframe with sample information
            includes columns: Name, group
        target_plot_df (pd.DataFrame): a dataframe with target sequence letters for plotting the left sequence heatmap - made by 'create_target_df_for_plotting'
            includes rows for each target, and columns for the target sequence letter to plot at that column (or .)
        output_folder (str): folder to write output to. If None, will be the current working directory
        file_prefix (str, optional): prefix for output folder Defaults to "CRISPRessoSea"
        crispresso2_info (dict, optional): crispresso2_info object for this CRISPRessoSea run (Defaults to None)
        heatmap_fig_height (int, optional): the height of the heatmap figure (Defaults to 24)
        heatmap_fig_width (int, optional): the width of the heatmap figure (Defaults to 24)
        heatmap_seq_plot_ratio (int, optional): the ratio of the width of the target sequence plot to the data plot (>1 means the seq plot is larger than the data plot) (Defaults to 1)
        heatmap_title_fontsize (int, optional): the fontsize of the plot titles for the heatmap (Defaults to 30)
        heatmap_y_tick_fontsize (int, optional): the fontsize of the y-axis tick labels for the heatmap (Defaults to 16)
        heatmap_x_tick_fontsize (int, optional): the fontsize of the x-axis tick labels for the heatmap (Defaults to 16)
        heatmap_nucleotide_fontsize (int, optional): the fontsize of the nucleotide labels in the heatmap(Defaults to 14)
        heatmap_legend_fontsize (int, optional): the fontsize of the legend title (Defaults to 20)
        heatmap_legend_ncol (int, optional): the number of columns in the legend (if None, each value legend value will have a columns)
        heatmap_max_value (float, optional): the maximum value for the heatmap color scale (Defaults to None)
        heatmap_min_value (float, optional): the minimum value for the heatmap color scale (Defaults to None)
        plot_group_order (list, optional): the order of the groups to plot (Defaults to None)
        dot_plot_ylims (list, optional): the y-axis limits [min,max] for the dot plot (Defaults to [None, None])
        sig_method_parameters (dict, optional): parameters for the significance method (Defaults to None)

    Returns:
        crispresso2_info with updated information about the added plot
    """
    gene_annotations = data_df["region_anno"].tolist()

    if output_folder is None:
        output_folder = os.getcwd()
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    groups = sample_df["group"].unique()
    groups = [group for group in groups if group != ""]  # remove empty groups
    groups = sorted(groups)

    plot_details = [
        {
            "col_to_plot": "highest_a_g_pct",
            "col_title": "Highest A->G %",
            "plot_suffix": "highest_a_g_pct",
        },
        {
            "col_to_plot": "highest_c_t_pct",
            "col_title": "Highest C->T %",
            "plot_suffix": "highest_c_t_pct",
        },
        {
            "col_to_plot": "highest_indel_pct",
            "col_title": "Highest Indel %",
            "plot_suffix": "highest_indel_pct",
        },
        {
            "col_to_plot": "tot_reads",
            "col_title": "Total Reads",
            "plot_suffix": "tot_reads",
        },
    ]

    for group in ["all", *groups]:
        cols_to_plot = [x["col_to_plot"] for x in plot_details]
        if group == "all":
            filtered_data_df = data_df
        else: # if plotting for a specific group, filter the data_df to only include the samples in that group
            sample_names = sample_df[sample_df.get("group", "") == group]["Name"].tolist()
            cols_to_keep = []
            for col in data_df.columns:
                if not any(col.endswith(plot_col) for plot_col in cols_to_plot):
                    cols_to_keep.append(col)
                else:
                    if any(col.startswith(sample_name) for sample_name in sample_names):
                        cols_to_keep.append(col)
            filtered_data_df = data_df[cols_to_keep]

        # for each type of plot (highest_a_g_pct, highest_c_t_pct, highest_indel_pct, tot_reads)
        for plot_detail in plot_details:
            col_to_plot = plot_detail["col_to_plot"]
            if group == "all":
                col_title = plot_detail["col_title"]
            else:
                col_title = f"{str(group).title()} {plot_detail['col_title']}"

            plot_suffix = f"{str(group)}_{plot_detail['plot_suffix']}"
            if file_prefix is not None:
                file_name = f"{file_prefix}_{plot_suffix}"
            else:
                file_name = plot_suffix

            plot_targets_and_heatmap(
                sample_df,
                target_plot_df,
                filtered_data_df,
                col_to_plot,
                col_title,
                row_annotations=gene_annotations,
                outfile_name=os.path.join(output_folder, file_name),
                fig_height=heatmap_fig_height,
                fig_width=heatmap_fig_width,
                seq_plot_ratio=heatmap_seq_plot_ratio,
                title_fontsize=heatmap_title_fontsize,
                y_tick_fontsize=heatmap_y_tick_fontsize,
                x_tick_fontsize=heatmap_x_tick_fontsize,
                nucleotide_fontsize=heatmap_nucleotide_fontsize,
                legend_title_fontsize=heatmap_legend_fontsize,
                df_data_heat_max=heatmap_max_value,
                df_data_heat_min=heatmap_min_value,
                legend_ncol=heatmap_legend_ncol,
                sig_method_parameters=sig_method_parameters,
            )
            if crispresso2_info is not None:
                crispresso2_info["results"]["general_plots"]["summary_plot_names"].append(file_name)
                crispresso2_info["results"]["general_plots"]["summary_plot_titles"][file_name] = col_title
                crispresso2_info["results"]["general_plots"]["summary_plot_labels"][file_name] = ""
                crispresso2_info["results"]["general_plots"]["summary_plot_datas"][file_name] = [("Aggregated Stats", "aggregated_stats_all.txt")]

            if group == "all":
                subset_data_df = data_df[ [ "target_id", "target_name", ] + [col for col in data_df.columns if col.endswith(col_to_plot)] ]
                dot_plot_suffix = f"{plot_suffix}_dot_plot"
                if file_prefix is not None:
                    file_name = f"{file_prefix}_{dot_plot_suffix}"
                else:
                    file_name = dot_plot_suffix

                plot_dot_plot(
                    subset_data_df,
                    col_to_plot,
                    col_title,
                    file_name,
                    sample_df,
                    output_folder,
                    plot_group_order=plot_group_order,
                    dot_plot_ylims=dot_plot_ylims,
                )
                if crispresso2_info is not None:
                    crispresso2_info["results"]["general_plots"][
                        "summary_plot_names"
                    ].append(file_name)
                    crispresso2_info["results"]["general_plots"]["summary_plot_titles"][
                        file_name
                    ] = (col_title + " Percent Edited")
                    crispresso2_info["results"]["general_plots"]["summary_plot_labels"][
                        file_name
                    ] = ""
                    crispresso2_info["results"]["general_plots"]["summary_plot_datas"][
                        file_name
                    ] = [("Aggregated Stats", "aggregated_stats_all.txt")]

    return crispresso2_info


# Function to melt DataFrame and plot dot plot
def plot_dot_plot(df, value_suffix, plot_title, file_prefix, sample_df, output_folder, plot_group_order=None, dot_plot_ylims=[None, None]):
    """
    Plots a dot plot showing average editing rate for each group

    Args:
        df (pd.DataFrame): DataFrame containing the columns 'target_id', 'target_name' and the columns to plot
        value_suffix (str): Suffix for the value column (e.g. 'highest_a_g_pct')
        plot_title (str): Title for the plot (e.g. 'Highest A->G %')
        file_prefix (str): Prefix for the output file name
        sample_df (pd.DataFrame): DataFrame containing sample information
        output_folder (str): Folder to save the output plot
        plot_group_order (list, optional): Order of groups to plot. If None, groups will be sorted alphabetically.
        dot_plot_ylims (list, optional): Y-axis limits [min,max] for the plot. Default is [None, None], which means auto-scaling.
    """
    group_counts = sample_df["group"].value_counts()
    
    melted_df = pd.melt(
        df,
        id_vars=[
            "target_id",
            "target_name",
        ],
        value_vars=[col for col in df.columns if col.endswith(value_suffix)],
        var_name="Sample",
        value_name="Value",
    )
    melted_df["Name"] = melted_df["Sample"].apply(
        lambda x: x.split(f"_{value_suffix}")[0]
    )
    melted_df.drop(columns=["Sample"], inplace=True)

    # melted_df now looks like this:
    # target_id    target_name  Name     Value
    # 0_Fs2_ON  guide1      sample1  0.5
    # 1_Fs2_OFF1  guide1      sample1  0.5

    if group_counts.shape[0] <= 1: # if only one group (or none)
        plt.figure(figsize=(20, 6))
        sns.stripplot(
            data=melted_df,
            x="target_name",
            y="Value",
            hue="Name",
            jitter=0.3,
            s=10,
            alpha=0.6,
        )
        plt.ylim(dot_plot_ylims[0], dot_plot_ylims[1])
        plt.title(f"Dot Plot of {plot_title} by Target")
        plt.xlabel("Target Name")
        plt.ylabel(plot_title)
        plt.xticks(fontsize=8, rotation=90)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        outfile_name = os.path.join(output_folder, file_prefix)
        plt.savefig(outfile_name + ".pdf", bbox_inches="tight")
        plt.savefig(outfile_name + ".png", bbox_inches="tight")
        plt.close()

    else:
        melted_df_with_group = melted_df.merge(sample_df, on='Name', how='left')
        """ this is some markdown"""

        """
        #has columns 
        # target_id   target_name   Value   Name    group
        # 0_f1_ON    gRNA_f1_ON   0.01    s1_ctl  ctl
        # 1_f1_OB1   gRNA_f1_OB1  0.01    s1_ctl  ctl
        """

        if plot_group_order is not None:
            melted_df_with_group['group'] = pd.Categorical(
                melted_df_with_group['group'], categories=plot_group_order, ordered=True
        )

        plt.figure(figsize=(20, 6))
        sns.stripplot(
                data=melted_df_with_group,
                x="target_name",
                y="Value",
                hue="group",
                jitter=0.3,
                s=5,
                alpha=0.6,
                palette="Dark2",
                dodge=True,
                legend=False
            )
        sns.barplot(
            data=melted_df_with_group,
            x="target_name",
            y="Value",
            hue="group",
            dodge=True,
            errorbar="sd",
            capsize=0.1,
            palette="Set2",
        )

        plt.ylim(dot_plot_ylims[0], dot_plot_ylims[1])
        plt.title(f"{plot_title} by Target and Group")
        plt.xlabel("Target Name")
        plt.ylabel(plot_title)
        plt.xticks(fontsize=8, rotation=90)
        plt.legend(title="Group")
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        outfile_name = os.path.join(output_folder, file_prefix)
        plt.savefig(outfile_name + ".pdf", bbox_inches="tight")
        plt.savefig(outfile_name + ".png", bbox_inches="tight")
        plt.close()


class SigMethod(Enum):
    NONE = 0
    HARD_CUTOFF = 1
    MEAN_DIFF = 2
    T_TEST = 3
    MANN_WHITNEY = 4
    NEG_BINOMIAL = 5

def parse_sig_method_parameters(sig_method_parameter_string, group_labels=None):
    """
    Parse the significance method parameters from a string or None.
    """
    significance_error_string = "Error in parsing sig_method_parameters value '"+sig_method_parameter_string+"': "
    significance_info_string = "Please provide a significance method description in the form of:\n\tnone\n\thard_cutoff,cutoff\n\tneg_binomial,group1,group2,alpha\n\tt_test,group1,group2,alpha\n\tmean_diff,group1,group2,cutoff"
    if sig_method_parameter_string is None or sig_method_parameter_string.strip() == "":
        return {
            'method': SigMethod.NONE,
        }
    
    parameter_els = sig_method_parameter_string.split(",")
    sig_method_parameters = {}
    if parameter_els[0].lower() == "none":
        return {
            'method': SigMethod.NONE,
        }
    elif parameter_els[0].lower() == "hard_cutoff" or parameter_els[0].lower() == "hardcutoff" or parameter_els[0].lower() == "hc":
        sig_method_parameters['method'] = SigMethod.HARD_CUTOFF
        if len(parameter_els) < 2:
            raise ValueError(significance_error_string + "Hard cutoff method requires a cutoff value.\n" + significance_info_string)
        try:
            sig_method_parameters['cutoff'] = float(parameter_els[1])
        except ValueError:
            raise ValueError(significance_error_string + "Cutoff value must be a number.")
    elif parameter_els[0].lower() == "mean_diff":
        sig_method_parameters['method'] = SigMethod.MEAN_DIFF
        if len(parameter_els) < 3:
            raise ValueError(significance_error_string + "Mean Difference method requires group_1, group_2, and cutoff values.\n" + significance_info_string)
        if group_labels is None or parameter_els[1] not in group_labels or parameter_els[2] not in group_labels:
            raise ValueError(significance_error_string + f"Groups {parameter_els[1]} and {parameter_els[2]} must be in the provided group labels: {group_labels}.")
        sig_method_parameters['group_1'] = parameter_els[1]
        sig_method_parameters['group_2'] = parameter_els[2]
        try:
            sig_method_parameters['cutoff'] = float(parameter_els[3]) if len(parameter_els) > 3 else 0.05
        except ValueError:
            raise ValueError(significance_error_string + "Cutoff value must be a number.")
    elif parameter_els[0].lower() == "mann_whitney" or parameter_els[0].lower() == "mw":
        sig_method_parameters['method'] = SigMethod.MANN_WHITNEY
        if len(parameter_els) < 4:
            raise ValueError(significance_error_string + "Mann-Whitney method requires group_1, group_2, and alpha values.\n" + significance_info_string)
        if group_labels is None or parameter_els[1] not in group_labels or parameter_els[2] not in group_labels:
            raise ValueError(significance_error_string + f"Groups {parameter_els[1]} and {parameter_els[2]} must be in the provided group labels: {group_labels}.")
        sig_method_parameters['group_1'] = parameter_els[1]
        sig_method_parameters['group_2'] = parameter_els[2]
        try:
            sig_method_parameters['alpha'] = float(parameter_els[3]) if len(parameter_els) > 3 else 0.05
        except ValueError:
            raise ValueError(significance_error_string + "Alpha value must be a number.")
    elif parameter_els[0].lower() == "t_test" or parameter_els[0].lower() == "ttest" or parameter_els[0].lower() == "t":
        sig_method_parameters['method'] = SigMethod.T_TEST
        if len(parameter_els) < 4:
            raise ValueError(significance_error_string + "T-test method requires group_1, group_2, and alpha values.\n" + significance_info_string)
        if group_labels is None or parameter_els[1] not in group_labels or parameter_els[2] not in group_labels:
            raise ValueError(significance_error_string + f"Groups {parameter_els[1]} and {parameter_els[2]} must be in the provided group labels: {group_labels}.")
        sig_method_parameters['group_1'] = parameter_els[1]
        sig_method_parameters['group_2'] = parameter_els[2]
        try:
            sig_method_parameters['alpha'] = float(parameter_els[3]) if len(parameter_els) > 3 else 0.05
        except ValueError:
            raise ValueError(significance_error_string + "Alpha value must be a number.")
    elif parameter_els[0].lower() == "neg_binomial" or parameter_els[0].lower() == "nbinom" or parameter_els[0].lower() == "nb":
        sig_method_parameters['method'] = SigMethod.NEG_BINOMIAL
        print('parameter_els', parameter_els)
        if len(parameter_els) < 4:
            raise ValueError(significance_error_string + "Negative Binomial method requires group_1, group_2, and alpha values.\n" + significance_info_string)
        if group_labels is None or parameter_els[1] not in group_labels or parameter_els[2] not in group_labels:
            raise ValueError(significance_error_string + f"Groups {parameter_els[1]} and {parameter_els[2]} must be in the provided group labels: {group_labels}.")
        sig_method_parameters['group_1'] = parameter_els[1]
        sig_method_parameters['group_2'] = parameter_els[2]
        try:
            sig_method_parameters['alpha'] = float(parameter_els[3]) if len(parameter_els) > 3 else 0.05
        except ValueError:
            raise ValueError(significance_error_string + "Alpha value must be a number.")
    else:
        raise ValueError(significance_error_string + f"Unknown significance method: {parameter_els[0]}. {significance_info_string}")
    
    return sig_method_parameters

def identify_significant_targets(df_data, df_total_count, sample_groups, sig_method_parameters):
    """
    Note that df_data values are percents - e.g. 0.5 = 0.5% editing, 25.5 = 25.5% editing
    Method 1: Hard cutoff, return values >= cutoff
    Method 2: Mean difference, return values where the mean difference between two groups is >= cutoff
    Method 3: Mann-Whitney, return values where the p-value is <= alpha
    Method 4: T-test, return values where the p-value is <= alpha
    Method 5: Negative Binomial, return values where the group effect is significant
     
    :param df_data: DataFrame with data to analyze, indexed by target_id and columns as samples
    :param df_total_count: DataFrame with total counts for each sample, indexed by sample name
    :param sample_groups: List of sample groups corresponding to the columns in df_data
    :param sig_method_parameters: Dictionary with parameters for the significance method, e.g. {'method': SigMethod.HARD_CUTOFF, 'cutoff': 0.05}

    :return: Tuple of DataFrame with True/False for significant results and a DataFrame with results table (if applicable)
    """
    
    if sig_method_parameters is None or sig_method_parameters['method'] == SigMethod.NONE:
        return None, None

    method = sig_method_parameters['method']
    if method == SigMethod.HARD_CUTOFF:
        cutoff = sig_method_parameters['cutoff']
        df_sig = df_data >= cutoff
        return df_sig, None

    elif method == SigMethod.MEAN_DIFF:
        df_sig = pd.DataFrame(False, index=df_data.index, columns=df_data.columns)

        res_table_rows = []
        group_1 = sig_method_parameters['group_1']
        group_1_inds = [i for i, x in enumerate(sample_groups) if x == group_1]
        group_2 = sig_method_parameters['group_2']
        group_2_inds = [i for i, x in enumerate(sample_groups) if x == group_2]

        # Don't do calculation if either group is empty
        if len(group_1_inds) == 0 or len(group_2_inds) == 0:
            return df_sig, None
        
        for idx, (df_idx, target_row) in enumerate(df_data.iterrows()):
            #check if any values are na
            if target_row[group_1_inds].isna().any() or target_row[group_2_inds].isna().any():
                res_table_rows.append([
                    group_1_data.nanmean(),
                    group_2_data.nanmean(),
                    np.nan,  # Coefficient for the group effect
                    'NA values present in data',
                    False,  # No significance if NA values are present
                ])
                continue
            group_1_data = target_row[group_1_inds].values
            group_2_data = target_row[group_2_inds].values

            # Calculate mean difference
            mean_diff = np.mean(group_1_data) - np.mean(group_2_data)
            is_greater_than_mean_diff = mean_diff >= sig_method_parameters['cutoff']

            # Append results to the list
            res_table_rows.append([
                group_1_data.mean(),
                group_2_data.mean(),
                mean_diff,
                'Mean Difference',
                is_greater_than_mean_diff
            ])
        res_table_df = pd.DataFrame(
            res_table_rows,
            columns=[group_1 + '_mean', group_2 + '_mean', 'mean_diff', 'test', 'significant']
        )
        res_table_df.index = df_data.index

        # Assign significant results to df_sig
        for idx, (df_idx, target_row) in enumerate(res_table_df.iterrows()):
            if target_row['significant']:
                if target_row[group_1 + '_mean'] > target_row[group_2 + '_mean']:
                    df_sig.iloc[idx, group_1_inds] = True
                else:
                    df_sig.iloc[idx, group_2_inds] = True

        #merge res_table_df with df_data
        res_table_df = res_table_df.merge(df_data, left_index=True, right_index=True)

        return df_sig, res_table_df

    elif method == SigMethod.MANN_WHITNEY:
        #create df of False with rows and columns of df_data
        df_sig = pd.DataFrame(False, index=df_data.index, columns=df_data.columns)

        res_table_rows = []
        group_1 = sig_method_parameters['group_1']
        group_1_inds = [i for i, x in enumerate(sample_groups) if x == group_1]
        group_2 = sig_method_parameters['group_2']
        group_2_inds = [i for i, x in enumerate(sample_groups) if x == group_2]

        #don't do calculation if either group is empty
        if len(group_1_inds) == 0 or len(group_2_inds) == 0:
            return df_sig, None

        for idx, target_row in df_data.iterrows():
            group_1_data = target_row[group_1_inds].values
            group_2_data = target_row[group_2_inds].values

            # Perform Mann-Whitney U test
            u_stat, p_val = sp.mannwhitneyu(group_1_data, group_2_data, alternative='two-sided')


            res_table_rows.append([
                group_1_data.mean(),
                group_2_data.mean(),
                u_stat,
                p_val,
                'Mann-Whitney U test',
            ])

        res_table_df = pd.DataFrame(res_table_rows, columns=[group_1 + '_mean', group_2 + '_mean', 'u_stat', 'p_val', 'test'])
        non_na_p_vals = res_table_df['p_val'].dropna()
        p_val_na_inds = res_table_df['p_val'].isna()

        p_corrected_res = multipletests(non_na_p_vals, method='fdr_bh', alpha=sig_method_parameters['alpha'])[1] #this is added to recent versions of scipy

        corrected_p = pd.Series(np.nan, index=res_table_df.index)
        corrected_p[~p_val_na_inds] = p_corrected_res

        res_table_df['bh_adj_p_val'] = corrected_p
        res_table_df['significant'] = res_table_df['bh_adj_p_val'] <= sig_method_parameters['alpha']
        res_table_df.index = df_data.index

        # .. hating how pandas does indexing and access with multi-indexes
        for idx, (df_idx, target_row) in enumerate(res_table_df.iterrows()):
            if target_row['significant']:
                if target_row[group_1 + '_mean'] > target_row[group_2 + '_mean']:
                    # Assign True to the appropriate columns for the current row
                    df_sig.iloc[idx,group_1_inds] = True
                else:
                    df_sig.iloc[idx,group_2_inds] = True

        res_table_df = res_table_df.merge(df_data, left_index=True, right_index=True)
        return df_sig, res_table_df

    elif method == SigMethod.T_TEST:
        #create df of False with rows and columns of df_data
        df_sig = pd.DataFrame(False, index=df_data.index, columns=df_data.columns)

        res_table_rows = []
        group_1 = sig_method_parameters['group_1']
        group_1_inds = [i for i, x in enumerate(sample_groups) if x == group_1]
        group_2 = sig_method_parameters['group_2']
        group_2_inds = [i for i, x in enumerate(sample_groups) if x == group_2]

        #don't do calculation if either group is empty
        if len(group_1_inds) == 0 or len(group_2_inds) == 0:
            return df_sig, None

        for idx, target_row in df_data.iterrows():
            group_1_data = target_row[group_1_inds].values
            group_2_data = target_row[group_2_inds].values

            t_stat, p_val = sp.ttest_ind(group_1_data, group_2_data, equal_var=False)

            res_table_rows.append([
                group_1_data.mean(),
                group_2_data.mean(),
                t_stat,
                p_val,
                'Welch T-test',
            ])

        res_table_df = pd.DataFrame(res_table_rows, columns=[group_1 + '_mean', group_2 + '_mean', 't_stat', 'p_val', 'test'])
        non_na_p_vals = res_table_df['p_val'].dropna()
        p_val_na_inds = res_table_df['p_val'].isna()

        p_corrected_res = multipletests(non_na_p_vals, method='fdr_bh', alpha=sig_method_parameters['alpha'])[1] #this is added to recent versions of scipy

        corrected_p = pd.Series(np.nan, index=res_table_df.index)
        corrected_p[~p_val_na_inds] = p_corrected_res

        res_table_df['bh_adj_p_val'] = corrected_p
        res_table_df['significant'] = res_table_df['bh_adj_p_val'] <= sig_method_parameters['alpha']
        res_table_df.index = df_data.index

        # .. hating how pandas does indexing and access with multi-indexes
        for idx, (df_idx, target_row) in enumerate(res_table_df.iterrows()):
            if target_row['significant']:
                if target_row[group_1 + '_mean'] > target_row[group_2 + '_mean']:
                    # Assign True to the appropriate columns for the current row
                    df_sig.iloc[idx,group_1_inds] = True
                else:
                    df_sig.iloc[idx,group_2_inds] = True

        res_table_df = res_table_df.merge(df_data, left_index=True, right_index=True)
        return df_sig, res_table_df




    elif method == SigMethod.NEG_BINOMIAL:
        """
        This is roughly equivalent to the R code:
        mod_counts = c(5,1,10,8,0,2,20,3)
        tot_reads = c(200, 950, 1200, 1100, 980, 1050, 1500, 1150)
        group = c(0,0,0,0,1,1,1,1)

        data = data.frame(tot_reads, mod_counts, group)
        data$mod_counts_adj = data$mod_counts + 1
        model <- glm.nb(mod_counts_adj ~ group + offset(log(tot_reads)), data = data, init.theta = 1)
        summary(model)

        ~~~~
        R output:
        ***
        group        -0.6510     0.5945  -1.095    0.274    
        ---
        Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

        (Dispersion parameter for Negative Binomial(1.8163) family taken to be 1)

            Null deviance: 9.7586  on 7  degrees of freedom
        Residual deviance: 8.6166  on 6  degrees of freedom
        AIC: 53.893

        Number of Fisher Scoring iterations: 1


                    Theta:  1.82 
                Std. Err.:  1.14 

        2 x log-likelihood:  -47.893

        ~~~

        Python output: 

        Generalized Linear Model Regression Results                  
        ==============================================================================
        Dep. Variable:         mod_counts_adj   No. Observations:                    8
        Model:                            GLM   Df Residuals:                        6
        Model Family:        NegativeBinomial   Df Model:                            1
        Link Function:                    Log   Scale:                          1.0000
        Method:                          IRLS   Log-Likelihood:                -24.394
        Date:                Mon, 02 Jun 2025   Deviance:                       5.4301
        Time:                        14:46:23   Pearson chi2:                     5.70
        No. Iterations:                     5   Pseudo R-squ. (CS):             0.1014
        Covariance Type:            nonrobust                                         
        ================================================================================
                        coef    std err          z      P>|z|      [0.025      0.975]
        --------------------------------------------------------------------------------
        const           -5.2140      0.539     -9.680      0.000      -6.270      -4.158
        group_binary     0.7193      0.761      0.945      0.345      -0.772       2.211
        ================================================================================
        
        R p-value: 0.274
        Python p-value: 0.345 
        
        #note that I can't get it to match up with Akcakaya et al. 2018, but it is close enough to be useful.
        #Note that this is the same as the R code:
        #mod_counts = c(17288, 13277, 13630, 625, 110, 91)
        #tot_reads = c(113499, 87435, 81144, 80495, 99696, 54993)
        #group = c(0,0,0,0,1,1)

        #supplemental 
        #rows 9-13
        mod_counts = c(17288, 13277, 13630, 625, 110, 91)
        tot_reads = c(113499, 87435, 81144, 80495, 99696, 54993)
        
        """
        # Create a DataFrame of False with the same shape as df_data
        df_sig = pd.DataFrame(False, index=df_data.index, columns=df_data.columns)

        res_table_rows = []
        group_1 = sig_method_parameters['group_1']
        group_1_inds = [i for i, x in enumerate(sample_groups) if x == group_1]
        group_2 = sig_method_parameters['group_2']
        group_2_inds = [i for i, x in enumerate(sample_groups) if x == group_2]

        # Don't do calculation if either group is empty
        if len(group_1_inds) == 0 or len(group_2_inds) == 0:
            return df_sig, None

        for idx, (df_idx, target_row) in enumerate(df_data.iterrows()):
            #check if any values are na
            group_1_data = target_row[group_1_inds].values / 100 # convert to proportions
            group_2_data = target_row[group_2_inds].values / 100

            if target_row[group_1_inds].isna().any() or target_row[group_2_inds].isna().any():
                res_table_rows.append([
                    np.nanmean(group_1_data)*100, #change back to percents
                    np.nanmean(group_2_data)*100,
                    np.nan,  # Coefficient for the group effect
                    np.nan,
                    'NA values present in data',
                ])
                continue

            group_1_total_counts = df_total_count.iloc[idx, group_1_inds].values
            group_2_total_counts = df_total_count.iloc[idx, group_2_inds].values
            group_1_mod_counts = group_1_data * group_1_total_counts
            group_2_mod_counts = group_2_data * group_2_total_counts

            group_labels = np.array([0] * len(group_1_data) + [1] * len(group_2_data))

            mod_counts = np.concatenate([group_1_mod_counts, group_2_mod_counts])
            total_counts = np.concatenate([group_1_total_counts, group_2_total_counts])

            df = pd.DataFrame({
                'group': group_labels,
                'total_counts': total_counts,
                'mod_counts': mod_counts
            })

            df['mod_counts_adj'] = df['mod_counts'] + 1

            # Design matrix
            df['group_binary'] = (df['group'] == 0).astype(int)
            X = sm.add_constant(df['group_binary'])

            # Response
            y = df['mod_counts_adj']

            # Offset (log total reads)
            offset = np.log(df['total_counts']+1)

            model = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=1), offset=offset)
            result = model.fit()

            p_val = result.pvalues['group_binary']

            # Combine data and create a group label
            combined_data = np.concatenate([group_1_data, group_2_data])
            combined_total_counts = np.concatenate([group_1_total_counts, group_2_total_counts])
            group_labels = np.array([0] * len(group_1_data) + [1] * len(group_2_data))

            res_table_rows.append([
                group_1_data.mean() * 100,  # Convert back to percentage
                group_2_data.mean() * 100,
                result.params['group_binary'],  # Coefficient for the group effect
                p_val,
                'Negative Binomial',
            ])

        res_table_df = pd.DataFrame(
            res_table_rows,
            columns=[group_1 + '_mean', group_2 + '_mean', 'coef', 'p_val', 'test']
        )
        non_na_p_vals = res_table_df['p_val'].dropna()
        p_val_na_inds = res_table_df['p_val'].isna()

        # Correct p-values using Benjamini-Hochberg FDR
        p_corrected_res = multipletests(non_na_p_vals, method='fdr_bh', alpha=sig_method_parameters['alpha'])[1]

        corrected_p = pd.Series(np.nan, index=res_table_df.index)
        corrected_p[~p_val_na_inds] = p_corrected_res

        res_table_df['bh_adj_p_val'] = corrected_p
        res_table_df['significant'] = res_table_df['bh_adj_p_val'] <= sig_method_parameters['alpha']
        res_table_df.index = df_data.index

        # Assign significant results to df_sig
        for idx, (df_idx, target_row) in enumerate(res_table_df.iterrows()):
            if target_row['significant']:
                if target_row[group_1 + '_mean'] > target_row[group_2 + '_mean']:
                    df_sig.iloc[idx, group_1_inds] = True
                else:
                    df_sig.iloc[idx, group_2_inds] = True

        res_table_df = res_table_df.merge(df_data, left_index=True, right_index=True)
        return df_sig, res_table_df

    return None, None


def plot_targets_and_heatmap(
    sample_df,
    target_plot_df,
    df_data,
    col_to_plot,
    df_data_title,
    df_data_heat_max=None,
    df_data_heat_min=None,
    outfile_name=None,
    row_annotations=None,
    fig_height=24,
    fig_width=24,
    target_names=None,
    seq_plot_ratio=1,
    title_fontsize=30,
    y_tick_fontsize=16,  # 16
    x_tick_fontsize=16,  # 16
    nucleotide_fontsize=14,  # 12
    legend_title_fontsize=20,  # 18
    legend_ncol=None,
    sig_method_parameters=None,
):
    """
    Plot a heatmap of target sequences (left) and the data (right)
    The plot of target sequences will contain
        '-' if there is a base missing in the target (relative to the on-target)
        '+' or 'AC' if there is a base added in the target (relative to the on-target)

    params:
    - sample_df: a dataframe of sample information with columns: Name, group
    - target_plot_df: a dataframe of target sequences with each column showing a nucleotide position, and each row showing the bases for each target
    - df_to_plot: a dataframe of data to plot in a heatmap with rows corresponding to the rows in target_plot_df
    - df_data_title: the title of the heatmap
    - df_data_heat_max: the maximum value for the heatmap
    - df_data_heat_min: the minimum value for the heatmap
    - outfile_name: the name of the output file (if None, the plot will be displayed)
    - row_annotations: a list of gene annotations for the rows in df_data (if None, no annotations will be displayed)
    - fig_height: the height of the figure
    - fig_width: the width of the figure
    - target_names: a list of target names to use for the y-axis of the heatmap (if None, the index of target_plot_df will be used)
    - seq_plot_ratio: the ratio of the width of the target sequence plot to the data plot (>1 means the seq plot is larger than the data plot)
    - title_fontsize: the fontsize of the plot titles
    - y_tick_fontsize: the fontsize of the y-axis tick labels
    - x_tick_fontsize: the fontsize of the x-axis labels
    - nucleotide_fontsize: the fontsize of the nucleotide labels
    - legend_title_fontsize: the fontsize of the legend title
    - legend_ncol: the number of columns in the legend (if None, each value legend value will have a columns)
    - sig_method_parameters: parameters for the significance method (if None, no significance will be calculated)

    returns:
    - None
    """
    cols_to_plot = [col for col in df_data.columns if col.endswith(col_to_plot)]
    if len(cols_to_plot) == 0:
        raise Exception("No columns found in df_data with suffix " + col_to_plot + " in columns: " + str(list(df_data.columns)))
    df_to_plot = df_data[cols_to_plot]
    df_to_plot.index = df_data['target_id']
    df_to_plot.columns = [col.replace("_" + col_to_plot, "") for col in cols_to_plot]

    #calculate total counts for negative binomial
    total_count_cols = [col for col in df_data.columns if col.endswith("tot_reads")]
    if len(total_count_cols) == 0:
        raise Exception("No total count columns found in " + str(list(df_data.columns)))
    df_total_count = df_data[total_count_cols]
    df_total_count.index = df_data['target_id']
    df_total_count.columns = [col.replace('_tot_reads', "") for col in cols_to_plot]

    sample_groups = []
    for sample_name in df_to_plot.columns:
        sample_group = sample_df[sample_df["Name"] == sample_name]["group"].values[0]
        sample_groups.append(sample_group)

    df_sig, df_pvals = identify_significant_targets(df_data=df_to_plot, df_total_count=df_total_count, sample_groups=sample_groups, sig_method_parameters=sig_method_parameters)

    if df_pvals is not None:
        p_value_output = outfile_name + ".p_values.txt"
        df_pvals.to_csv(p_value_output, sep="\t", index=True, header=True)
        info('Saved p-values to ' + p_value_output)

    # Create a custom color palette for the letters
    color_mapping = {
        "A": "#E3EFA9",
        "T": "#CCCFE0",
        "C": "#FBC6C6",
        "G": "#FCE588",
        ".": "#F4F4F6",
        "-": "#BEC0C6",
        "+": "#BEC0C6",
        "default": "gray",
    }
    # color_mapping = {'A': '#90adc6', 'T': '#e9eaec', 'C': '#fad02c', 'G': '#333652', '.': '#c8df52', '-':'#F67E7D',  'default': 'gray'}
    # color_mapping = {'A': '#FFB7B2', 'T': '#C3CDE6', 'C': '#E2F0CB', 'G': '#B5EAD7', '.': '#FFDAC1', 'default': 'gray'}
    unique_letters = color_mapping.keys()
    colors = [color_mapping[letter] for letter in unique_letters]
    cmap = ListedColormap(colors)

    # Map letters to integers for heatmap
    letter_to_int = {letter: i for i, letter in enumerate(unique_letters)}

    def get_int_from_letter(letter):
        if letter in letter_to_int:
            return letter_to_int[letter]
        if len(letter) > 1:
            return letter_to_int["+"]
        return letter_to_int["default"]

    df_mapped = target_plot_df.applymap(lambda x: get_int_from_letter(x))

    # Create the figure and gridspec
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[10, 1],
        width_ratios=[seq_plot_ratio, 1],
        hspace=0,
        wspace=0.05,
    )

    # Create the first heatmap (target seqs and mismatches)
    ax1 = fig.add_subplot(gs[0, 0])
    sns.heatmap(
        df_mapped,
        annot=target_plot_df,
        fmt="",
        cmap=cmap,
        cbar=False,
        ax=ax1,
        xticklabels=False,
        vmin=0,
        vmax=len(unique_letters) - 1,
        annot_kws={"fontsize": nucleotide_fontsize},
    )
    ax1.tick_params(axis="y", labelsize=y_tick_fontsize)
    # replace y tick lables with target names
    if target_names is not None:
        ax1.set_yticklabels(target_names, rotation=0, fontsize=y_tick_fontsize)
    ax1.set_title("Target sequences", fontsize=title_fontsize)
    ax1.set_ylabel('')

    # Create the second heatmap (continuous data from df_data)
    ax2 = fig.add_subplot(gs[0, 1])
    data_annot = True
    if df_to_plot.shape[0] > 20:
        data_annot = False
    if df_data_heat_max is None:
        df_data_heat_max = df_to_plot.max(skipna=True).max(skipna=True)
    if df_data_heat_min is None:
        df_data_heat_min = df_to_plot.min(skipna=True).min(skipna=True)

    if row_annotations is not None:
        row_annotations = ["" if pd.isna(anno) else anno for anno in row_annotations]

    if row_annotations is not None and set(row_annotations) != {
        ""
    }:  # adjust padding if no annotations
        max_anno = max([len(anno) for anno in row_annotations]) / 100
        max_anno = max(max_anno, 0.03)
        cbar_pad = max_anno
        # cbar_pad = 0.11
    else:
        cbar_pad = 0.03
    ax_heat = sns.heatmap(
        df_to_plot,
        annot=data_annot,
        cmap="Blues",
        ax=ax2,
        yticklabels=False,
        vmin=df_data_heat_min,
        vmax=df_data_heat_max,
        cbar_kws={"orientation": "vertical", "pad": cbar_pad},
    )
    ax_heat.collections[0].cmap.set_bad("0.7")
    ax_heat.collections[0].colorbar.ax.tick_params(labelsize=y_tick_fontsize)
    ax_heat.set_ylabel('')
    ax2.set_title(df_data_title, fontsize=title_fontsize)
    ax2.tick_params(axis="x", labelsize=x_tick_fontsize, labelrotation=90)
    ax2.set_yticks([])
    ax2.set_ylabel('')

    info("Plotted heatmap to " + outfile_name)

    # Adding gene annotation labels to right of heatmap
    if row_annotations is not None and set(row_annotations) != {
        ""
    }:  # if all annotations are empty, don't plot them
        ax2.set_yticks(
            np.arange(len(row_annotations)) + 0.5,
            labels=row_annotations,
            fontsize=y_tick_fontsize,
        )
        ax2.tick_params(
            axis="y",
            which="both",
            left=False,
            right=True,
            labelright=True,
            labelleft=False,
        )

    # Annotation of significant values (red box)
    if df_sig is not None:
        for i in range(df_sig.shape[0]):
            for j in range(df_sig.shape[1]):
                if df_sig.iloc[i, j]:
                    ax_heat.add_patch(
                        plt.Rectangle(
                            (j, i), 1, 1, fill=False, edgecolor="red", lw=2, alpha=0.6
                        )
                    )

    # Create a custom legend under the first heatmap
    ax_legend = fig.add_subplot(gs[1, 0])
    ax_legend.axis("off")
    patches = [
        plt.plot(
            [],
            [],
            marker="s",
            ms=10,
            ls="",
            mec=None,
            color=color_mapping[letter],
            label="{:s}".format(letter),
        )[0]
        for letter in color_mapping
    ]

    if legend_ncol is not None:
        legend_ncol = legend_ncol
    else:
        legend_ncol = len(color_mapping)

    ax_legend.legend(
        handles=patches,
        loc="center",
        title="Target nucleotides",
        ncol=legend_ncol,
        bbox_to_anchor=(0.5, 1),
        title_fontsize=legend_title_fontsize,
        fontsize=legend_title_fontsize - 2,
    )

    gs.tight_layout(fig)

    if outfile_name is not None:
        plt.savefig(outfile_name + ".pdf", bbox_inches="tight")
        plt.savefig(outfile_name + ".png", bbox_inches="tight")
    else:
        plt.show()

    plt.close(fig)


def add_region_annotations_to_target_df(target_df, gene_annotations):
    """
    Add gene annotations to the target dataframe

    params:
    - target_df: the target dataframe. This function expects the columns: target_chr, target_pos
    - gene_annotations: the path to the gene annotations file. This function expects a column for:
        chromosome ('chrom' or 'chr'),
        start ('txStart','start', or 'Start'),
        end ('txEnd', 'end', or 'End'), and
        gene name ('name' and/or 'name2')

    returns:
    - the target dataframe with gene annotations added in the column 'region_anno'
    """
    info(f"Loading gene coordinates from annotation file: {gene_annotations}...")
    try:
        df_genes = pd.read_csv(gene_annotations, sep="\t")

        if "chr" in df_genes.columns:
            df_genes["chrom"] = df_genes.chr
        if "txStart" in df_genes.columns:
            df_genes["txStart"] = df_genes.txStart.astype(int)
        elif "start" in df_genes.columns:
            df_genes["txStart"] = df_genes.start.astype(int)
        elif "Start" in df_genes.columns:
            df_genes["txStart"] = df_genes.Start.astype(int)
        else:
            raise Exception(
                'Cannot find start column in gene annotations file. Expecting columns "txStart" or "start". Only observed columns '
                + str(df_genes.columns)
            )

        if "txEnd" in df_genes.columns:
            df_genes["txEnd"] = df_genes.txEnd.astype(int)
        elif "end" in df_genes.columns:
            df_genes["txEnd"] = df_genes.end.astype(int)
        elif "End" in df_genes.columns:
            df_genes["txEnd"] = df_genes.End.astype(int)
        else:
            raise Exception(
                'Cannot find end column in gene annotations file. Expecting columns "txEnd" or "end". Only observed columns '
                + str(df_genes.columns)
            )

        df_genes.head()
    except Exception:
        df_genes = pd.DataFrame(
            columns=[
                "chrom",
                "txStart",
                "txEnd",
                "name",
            ]
        )
        warn("Failed to load the gene annotations file.")

    def find_overlapping_genes(row, df_genes):
        df_genes_overlapping = df_genes.loc[
            (df_genes.chrom == row.target_chr)
            & (df_genes.txStart <= row.target_pos)
            & (row.target_pos <= df_genes.txEnd)
        ]
        genes_overlapping = []

        for idx_g, row_g in df_genes_overlapping.iterrows():
            if "name" in row_g.keys() and "name2" in row_g.keys():
                genes_overlapping.append("%s (%s)" % (row_g.name2, row_g["name"]))
            elif "#name" in row_g.keys() and "name2" in row_g.keys():
                genes_overlapping.append("%s (%s)" % (row_g.name2, row_g["#name"]))
            elif "#name" in row_g.keys():
                genes_overlapping.append("%s" % (row_g["#name"]))
            elif "name" in row_g.keys():
                genes_overlapping.append("%s" % (row_g["name"]))
            else:
                genes_overlapping.append("%s" % (row_g[0]))

        row["region_anno"] = ",".join(genes_overlapping)
        return row

    target_df["region_anno"] = ""
    target_df = target_df.apply(lambda row: find_overlapping_genes(row, df_genes), axis=1)

    return target_df


def replot(
    reordered_stats_file,
    reordered_sample_file,
    output_folder=None,
    file_prefix=None,
    name_column=None,
    fig_height=24,
    fig_width=24,
    seq_plot_ratio=1,
    title_fontsize=30,
    y_tick_fontsize=16,  # 16
    x_tick_fontsize=16,  # 16
    nucleotide_fontsize=14,  # 12
    legend_title_fontsize=20,  # 18
    legend_ncol=None,
    plot_group_order=None,
    dot_plot_ylims=[None, None],
    heatmap_max_value=None,
    heatmap_min_value=None,
    sig_method_parameters=None,
):
    """
    Replot a completed analysis using a reordered target file

    params:
    - reordered_stats_file: the reordered stats file based on (having the same columns and layout as) a previously-completed aggregated_stats_all.txt file
    - reordered_sample_file: path to the sample file with headers: Name, group, fastq_r1, fastq_r2 (group is always optional, fastq_r2 is optional for single-end reads)
    - output_folder: the output folder for the output plots
    - file_prefix: the prefix for the output plots
    - name_column: the column to use as the displayed name for each sample in the plot
    - fig_height: the height of the figure
    - fig_width: the width of the figure
    - seq_plot_ratio: the ratio of the width of the target sequence plot to the data plot (>1 means the seq plot is larger than the data plot)
    - title_fontsize: the fontsize of the plot titles for the heatmap
    - y_tick_fontsize: the fontsize of the y-axis tick labels for the heatmap
    - x_tick_fontsize: the fontsize of the x-axis tick labels for the heatmap
    - nucleotide_fontsize: the fontsize of the nucleotide labels for the heatmap
    - legend_title_fontsize: the fontsize of the legend title for the heatmap
    - legend_ncol: the number of columns in the legend for the heatmap (if None, each value legend value will have a columns)
    - plot_group_order: the order of the groups to plot (if None, the groups will be sorted alphabetically)
    - heatmap_max_value: the maximum value for the heatmap color scale (if None, the maximum value will be determined automatically)
    - heatmap_min_value: the minimum value for the heatmap color scale (if None, the minimum value will be determined automatically)
    - sig_method_parameters: parameters for the significance method (if None, no significance will be calculated)
    """

    sample_df = parse_sample_file(reordered_sample_file)

    if output_folder is not None:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder, exist_ok=True)

    reordered_stats_df = pd.read_csv(reordered_stats_file, sep="\t")
    target_plot_df = create_target_df_for_plotting(reordered_stats_df)

    if name_column is not None:
        if name_column not in sample_df.columns:
            raise Exception(f"Column '{name_column}' not found in sample file. Available columns: {list(sample_df.columns)}")
            

    info("Plotting for " + str(len(reordered_stats_df)) + " targets", {"percent_complete": 20})

    # reorder the stats column to the order of the sample file
    sample_cols_reordered = []
    for sample in sample_df["Name"]:
        for col in reordered_stats_df.columns:
            if col.startswith(sample):
                sample_cols_reordered.append(col)

    reordered_cols = []
    for col in reordered_stats_df.columns:
        if col not in sample_cols_reordered:
            reordered_cols.append(col)

    reordered_cols += sample_cols_reordered
    reordered_stats_df = reordered_stats_df[reordered_cols]

    create_plots(
        data_df=reordered_stats_df,
        sample_df=sample_df,
        target_plot_df=target_plot_df,
        output_folder=output_folder,
        file_prefix=file_prefix,
        heatmap_fig_height=fig_height,
        heatmap_fig_width=fig_width,
        heatmap_seq_plot_ratio=seq_plot_ratio,
        heatmap_title_fontsize=title_fontsize,
        heatmap_y_tick_fontsize=y_tick_fontsize,
        heatmap_x_tick_fontsize=x_tick_fontsize,
        heatmap_nucleotide_fontsize=nucleotide_fontsize,
        heatmap_legend_fontsize=legend_title_fontsize,
        heatmap_legend_ncol=legend_ncol,
        heatmap_max_value=heatmap_max_value,
        heatmap_min_value=heatmap_min_value,
        plot_group_order=plot_group_order,
        dot_plot_ylims=dot_plot_ylims,
        sig_method_parameters=sig_method_parameters,
    )

    info('Replotting complete', {"percent_complete": 100})

def make_target_info_file(guide_seq_str, guide_name_str, guide_pam, genome_file, max_mismatches, max_rna_bulges, max_dna_bulges, output_folder, file_prefix):
    """
    Make a guide info file from a list of guides. The guide info file will contain the columns:
    If not assigned, guide names will be assigned as "guide_0", "guide_1", etc.

    params:
    - guide_seq_str: a string of comma-separated guide sequences
    - guide_name_str: a string of comma-separated guide names
    - guide_pam: the PAM sequence
    - genome_file: the genome file for input to casoffinder
    - max_mismatches: the maximum number of mismatches between guide and off-target to search for
    - max_rna_bulges: the maximum number of RNA bulges between guide and off-target to search for (Note that Cas-OFFinder 3 is required to detect bulges in off-targets)
    - max_dna_bulges: the maximum number of DNA bulges between guide and off-target to search for (Note that Cas-OFFinder 3 is required to detect bulges in off-targets)
    - output_folder: the output folder
    - file_prefix: the prefix for the output file

    """

    info('Creating guide info file', {"percent_complete": 0})

    casoffinder_input_file = output_folder + "/" + file_prefix + ".casoffinder_input.txt"
    casoffinder_output_file = output_folder + "/" + file_prefix + ".casoffinder_output.txt"
    casoffinder_log_file = output_folder + "/" + file_prefix + ".casoffinder_log.txt"
    casoffinder_annotations = []

    guide_output_file = output_folder + "/" + file_prefix + ".guide_info.txt"


    #link in genome -- it throws a seg fault if it's in the bowtie directory
    linked_genome = os.path.abspath(output_folder + '.genome.fa')
    if os.path.exists(linked_genome):
        os.remove(linked_genome)
    os.symlink(os.path.abspath(genome_file),linked_genome)

    valid_nucs = ['A','C','G','T','N']
    guide_seqs = guide_seq_str.split(",")
    guide_names = guide_name_str.split(",")

    for guide_seq in guide_seqs:
        for nuc in guide_seq:
            if nuc.upper() not in valid_nucs:
                raise Exception('Invalid nucleotide found in guide sequence "%s": %s'%(guide_seq,nuc))

    while len(guide_seqs) > len(guide_names):
        curr_guide_id = 0
        potential_guide_name = "guide_%d" % curr_guide_id
        if potential_guide_name not in guide_names:
            guide_names.append(potential_guide_name)
        curr_guide_id += 1

    guide_len = max([len(x) for x in guide_seqs])
    pam_len = len(guide_pam)
    if not os.path.exists(casoffinder_input_file):
        with open (casoffinder_input_file,'w') as cout:
            cout.write(linked_genome+"\n")
            if max_dna_bulges > 0 or max_rna_bulges > 0: ## max bulges should be written for Cas-OFFinder 3 but they mess up the output for Cas-OFFinder 2
                cout.write("N"*guide_len + guide_pam + " " + str(max_dna_bulges) + " " + str(max_rna_bulges) + "\n")
            else:
                cout.write("N"*guide_len + guide_pam + "\n")
            for guide in guide_seqs:
                cout.write(guide + "N"*pam_len + " " + str(max_mismatches) + "\n")

    casoffinder_cmd = '(%s %s C %s) &>> %s'%('cas-offinder',casoffinder_input_file,casoffinder_output_file,casoffinder_log_file)
    ## tell the user to run this simplified command (without redirecting output to log file) if we can't run it here
    external_casoffinder_cmd = '%s %s C %s'%('cas-offinder',casoffinder_input_file,casoffinder_output_file)

    debug('Creating casoffinder command: "%s"'%casoffinder_cmd, {"percent_complete": 20})

    with open (casoffinder_log_file,'w') as cout:
        cout.write('Linking genome from %s to %s\n'%(genome_file,linked_genome))
        cout.write('Command used:\n===\n%s\n===\nOutput:\n===\n'%casoffinder_cmd)

    debug(casoffinder_cmd)
    if not os.path.exists(casoffinder_output_file):
        debug(f'Cannot find casoffinder output file: {casoffinder_output_file}.')
        #check casoffinder
        try:
            casoffinder_result = subprocess.check_output('cas-offinder --version', stderr=subprocess.STDOUT, shell=True)
        except Exception:
            raise Exception('Error: Cas-OFFinder is required. Please run the command: "' + external_casoffinder_cmd + '" manually and then rerun')
        if not 'Cas-OFFinder v' in str(casoffinder_result):
            raise Exception('Error: Cas-OFFinder is required. Please run the command: "' + external_casoffinder_cmd + '" manually and then rerun')

        casoffinder_result = subprocess.check_output(casoffinder_cmd,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        debug('Casoffinder output:' + casoffinder_result)
    else:
        debug('Using previously-calculated offtargets')

    os.remove(linked_genome)
    info('Parsing casoffinder output', {"percent_complete": 60})

    is_casoffinder_3 = False
    with open(casoffinder_output_file,'r') as cin:
        line1 = cin.readline()
        if line1.startswith('##Generated by Cas-OFFinder 3'):
            is_casoffinder_3 = True
    if (max_dna_bulges > 0 or max_rna_bulges > 0) and not is_casoffinder_3:
        warn('Note that Cas-OFFinder 3 is required to detect bulges in off-targets. Cas-OFFinder output is detected to be from Cas-OFFinder version 2. Proceeding with only mismatches.')

    if is_casoffinder_3:
        with open(casoffinder_output_file,'r') as cin:
            info_line = cin.readline() # ##Generated by Cas-OFFinder
            head_line = cin.readline()
            head_line_els = head_line.strip().split("\t")
            if head_line_els[1] != 'Bulge Type':
                raise Exception('Cannot parse casoffinder output. Expecting "Bulge Type" in second column of header. Found: ' + head_line)
            if head_line_els[4] != 'Chromosome':
                raise Exception('Cannot parse casoffinder output. Expecting "Chromosome" in fifth column of header. Found: ' + head_line)

            for line in cin:
                line_els = line.strip().split("\t")
                ontarget_seq = line_els[2][0:-pam_len]
                chrom = line_els[4]
                loc = line_els[5]
                strand = line_els[6]
                mismatch_count = int(line_els[7])
                bulge_count = int(line_els[8])
                total_diff_count = mismatch_count + bulge_count

                target_name = None
                for guide_seq in guide_seqs:
                    if ontarget_seq.replace("-","").replace("N","").upper() == guide_seq.upper():
                        target_name = guide_names[guide_seqs.index(guide_seq)]
                        break
                if target_name is None:
                    warn('Cannot parse casoffinder output. Cannot find guide name for guide sequence: ' + ontarget_seq)
                    target_name = ontarget_seq.replace("-","").replace("N","")
                
                offtarget_seq = line_els[3][1:-pam_len].replace("-","").replace("N","")
                offtarget_pam = line_els[3][-pam_len:]

                mismatch_string = 'Mismatches: ' + str(mismatch_count) + ', Bulges: ' + line_els[1] + ' ' + str(bulge_count)

                locus_string = chrom+":"+strand+loc
                casoffinder_annotations.append([target_name, target_name + "_OB" + str(total_diff_count) + "_" + chrom + "_" + loc,
                                                     offtarget_seq, offtarget_pam, str(total_diff_count), locus_string, mismatch_string])
    else: # parse old casoffinder output
        with open (casoffinder_output_file,'r') as cin:
            for line in cin:
                line_els = line.strip().split("\t")
                seq = line_els[0][0:-pam_len]
                chrom = line_els[1]
                loc = line_els[2]
                strand = line_els[4]
                mismatch_count = line_els[5]
                total_diff_count = mismatch_count

                offtarget_seq = line_els[3][0:-pam_len]
                offtarget_pam = line_els[3][-pam_len:]

                target_name = None
                for guide_seq in guide_seqs:
                    if seq.replace("-","").upper() == guide_seq.upper():
                        target_name = guide_names[guide_seqs.index(guide_seq)]
                        break
                if target_name is None:
                    warn('Cannot parse casoffinder output. Cannot find guide name for guide sequence: ' + seq)
                    target_name = guide_seq

                mismatch_string = 'Mismatches: ' + str(mismatch_count)

                locus_string = chrom+":"+strand+loc
                casoffinder_annotations.append([target_name, target_name + "_OB" + str(total_diff_count) + "_" + chrom + "_" + loc,offtarget_seq, offtarget_pam, str(total_diff_count), locus_string, mismatch_string])

    parsed_output = pd.DataFrame(casoffinder_annotations,columns=['Guide','Target','Sequence','PAM','#MM','Locus','Mismatch_info'])
    info('Found %d off-target sites using Cas-OFFinder'%(parsed_output.shape[0]), {'percent_complete': 80})

    parsed_output = parsed_output.sort_values(by=['#MM'], ascending=True)

    parsed_output.to_csv(guide_output_file, sep="\t", index=False)

    info('Wrote guide info to ' + guide_output_file, {'percent_complete': 100})

# main entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process multiple pooled sequencing runs",
    )
    
    subparsers = parser.add_subparsers(
        dest="subcommand", help="Process a new run or Replot a completed run"
    )
    ### Process
    process_parser = subparsers.add_parser(
        "Process",
        help="Process a new set of samples",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    process_parser.add_argument(
        "-o", "--output_folder", help="Output folder", default=None
    )
    process_parser.add_argument(
        "-t",
        "--target_file",
        help="Target file - list of targets with one target per line and the headers ['Guide','Sequence','PAM','#MM','Locus']",
        required=True,
    )
    process_parser.add_argument(
        "-s",
        "--sample_file",
        help="Sample file - list of samples with one sample per line with headers ['Name','fastq_r1','fastq_r2']",
        required=True,
    )
    process_parser.add_argument(
        "--gene_annotations",
        help='Gene annotations file - a tab-separated .bed file with gene annotations with columns "chr" or "chrom", "start" or "txstart", and "end" or "txend" as well as "name"',
        default=None,
    )
    process_parser.add_argument(
        "-x",
        "--genome_file",
        help="Bowtie2-indexed genome file - files ending in and .bt2 must be present in the same folder.",
        required=True,
    )
    process_parser.add_argument(
        "-p",
        "--n_processes",
        help='Number of processes to use. Set to "max" to use all available processors.',
        type=str,
        default="8",
    )
    process_parser.add_argument(
        "--crispresso_quantification_window_center",
        help="Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17 to only include mutations near the 5' end of the sgRNA.",
        default=-3,
        type=int,
    )
    process_parser.add_argument(
        "--crispresso_quantification_window_size",
        help="Size (in bp) of the quantification window extending from the position specified by the '--cleavage_offset' or '--quantification_window_center' parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp.",
        default=1,
        type=int,
    )
    process_parser.add_argument(
        "--crispresso_base_editor_output",
        help='Outputs plots and tables to aid in analysis of base editor studies.',
        action="store_true",
    )
    process_parser.add_argument(
        "--crispresso_default_min_aln_score",
        help="Default minimum homology score for a read to align to a reference amplicon.",
        default=60,
        type=int,
    )
    process_parser.add_argument(
        "--crispresso_plot_window_size",
        help="Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the quantification_window_center for each guide are plotted.",
        default=20,
        type=int,
    )

    process_parser.add_argument(
        "--allow_unplaced_chrs",
        help="Allow regions on unplaced chromosomes (chrUn, random, etc). By default, regions on these chromosomes are excluded. If set, regions on these chromosomes will be included.",
        action="store_true",
    )
    process_parser.add_argument(
        "--plot_only_complete_targets",
        help="Plot only targets with all values. If not set, all targets will be plotted.",
        action="store_true",
    )
    process_parser.add_argument(
        "--min_amplicon_coverage",
        help="Minimum number of reads to cover a location for it to be plotted. Otherwise, it will be set as NA",
        default=10,
        type=int,
    )
    process_parser.add_argument(
        "--sort_based_on_mismatch",
        help="Sort targets based on mismatch count. If true, the on-target will always be first",
        action="store_true",
    )
    process_parser.add_argument(
        "--allow_target_match_to_other_region_loc",
        help="If true, targets can match to regions even if the target chr:start is not in that region (e.g. if the target sequence is found in that region). If false/unset, targets can only match to regions matching the target chr:start position. This flag should be set if the genome for guide design was not the same as the analysis genome.",
        action="store_true",
    )
    process_parser.add_argument(
        "--top_percent_cutoff",
        help="The top percent of aligned regions (by region read depth) to consider in finding non-overlapping regions during demultiplexing. This is a float between 0 and 1. For example, if set to 0.2, the top 20%% of regions (by read depth) will be considered.",
        default=0.2,
        type=float,
    )
    process_parser.add_argument(
        "--min_amplicon_len",
        help="The minimum length of an amplicon to consider in finding non-overlapping regions during demultiplexing. Amplicons shorter than this will be ignored.",
        type=int,
        default=50,
    )
    process_parser.add_argument(
        "--fail_on_pooled_fail",
        help="If true, fail if any pooled CRISPResso run fails. By default, processing will continue even if sub-CRISPResso commands fail.",
        action="store_true",
    )
    process_parser.add_argument(
        "--plot_group_order",
        help="Order of the groups to plot (if None, the groups will be sorted alphabetically)",
        default=None,
        type=str,
    )
    process_parser.add_argument(
        "--sig_method_parameters",
        help="Parameters for the significance method in the form of:\n\tnone\n\thard_cutoff,cutoff\n\tneg_binomial,group1,group2,alpha\n\tt_test,group1,group2,alpha\n\tmean_diff,group1,group2,cutoff",
        default=None,
        type=str,
    )
    process_parser.add_argument(
        "-v",
        "--verbosity",
        help="Verbosity level of output to the console (1-4) 4 is the most verbose",
        type=int,
        default=3,
    )
    process_parser.add_argument(
        "--debug",
        help="Print debug information",
        action="store_true",
    )

    ### Replot
    plot_parser = subparsers.add_parser(
        "Replot",
        help="Replot completed analysis using reordered sample table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plot_parser.add_argument(
        "-o", "--output_folder", help="Output folder where output plots should be placed", default=None
    )
    plot_parser.add_argument(
        "-p",
        "--file_prefix",
        help="File prefix for output files",
        default="CRISPRessoSea",
    )
    plot_parser.add_argument(
        "-f",
        "--reordered_stats_file",
        help="Reordered statistics file - made by reordering rows from aggregated_stats_all.txt",
        required=True,
    )
    plot_parser.add_argument(
        "-s",
        "--reordered_sample_file",
        help="Reordered_sample_file - path to the sample file with headers: Name, group, fastq_r1, fastq_r2 (group is always optional, fastq_r2 is optional for single-end reads)",
        required=True,
    )
    plot_parser.add_argument(
        "-n",
        "--name_column",
        help="Column name to set as the displayed name for each sample in the plot",
        default=None,
    )
    plot_parser.add_argument(
        "--fig_width", help="Width of the figure", default=24, type=int
    )
    plot_parser.add_argument(
        "--fig_height", help="Height of the figure", default=24, type=int
    )
    plot_parser.add_argument(
        "--seq_plot_ratio",
        help="Ratio of the width of the sequence plot to the data plot (>1 means the seq plot is larger than the data plot)",
        default=1,
        type=float,
    )
    plot_parser.add_argument(
        "--plot_group_order",
        help="Order of the groups to plot (if None, the groups will be sorted alphabetically)",
        default=None,
        type=str,
    )
    plot_parser.add_argument(
        "--sig_method_parameters",
        help="Parameters for the significance method in the form of:\n\tnone\n\thard_cutoff,cutoff\n\tneg_binomial,group1,group2,alpha\n\tt_test,group1,group2,alpha\n\tmean_diff,group1,group2,cutoff",
        default=None,
        type=str,
    )
    plot_parser.add_argument(
        "--dot_plot_ylims",
        help="Comma-separated min,max y-axis limits for the dot plot. If None, the y-axis limits will be set automatically.",
        default='None,None',
        type=str,
    )
    plot_parser.add_argument(
        "--heatmap_max_value",
        help="Maximum value for the heatmap color scale, where a value of 1 sets the max value color to 1% (if None, the maximum value will be determined automatically)",
        default=None,
        type=float,
    )
    plot_parser.add_argument(
        "--heatmap_min_value",
        help="Minimum value for the heatmap color scale, where a value of 1 sets the min value color to 1% (if None, the minimum value will be determined automatically)",
        default=None,
        type=float,
    )
    plot_parser.add_argument(
        "-v",
        "--verbosity",
        help="Verbosity level of output to the console (1-4) 4 is the most verbose",
        type=int,
        default=3,
    )
    plot_parser.add_argument(
        "--debug",
        help="Print debug information",
        action="store_true",
    )

    ### MakeGuideFile
    makeguidefile_parser = subparsers.add_parser(
        "MakeGuideFile",
        help="Make a guide file by computationally enumerating off-targets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    makeguidefile_parser.add_argument(
        "-o", "--output_folder", help="Output folder", default=None
    )
    makeguidefile_parser.add_argument(
        "-p",
        "--file_prefix",
        help="File prefix for output files",
        default="CRISPRessoSea",
    )
    makeguidefile_parser.add_argument(
        "-g",
        "--guide_seq",
        help="Guide sequence(s) to create a guide file for. Multiple guides may be separated by commas.",
        required=True,
    )
    makeguidefile_parser.add_argument(
        "-gn",
        "--guide_name",
        help="Guide name(s) to use for the guide sequence(s). Multiple names may be separated by commas.",
        default='guide_0'
    )
    makeguidefile_parser.add_argument(
        "--pam",
        help="PAM sequence(s) to use for the guide sequence(s).",
        default='NGG'
    )
    makeguidefile_parser.add_argument(
        "-x",
        "--genome_file",
        help="Bowtie2-indexed genome file - files ending in and .bt2 must be present in the same folder.",
        required=True,
    )
    makeguidefile_parser.add_argument(
        "--max_mismatches",
        help="Maximum number of mismatches to allow in the discovered offtargets",
        default=4,
        type=int,
    )
    makeguidefile_parser.add_argument(
        "--max_dna_bulges",
        help="Maximum number of DNA bulges to allow in the discovered offtargets. Note that Cas-OFFinder 3 is required to detect sites with bulges.",
        default=0,
        type=int,
    )
    makeguidefile_parser.add_argument(
        "--max_rna_bulges",
        help="Maximum number of RNA bulges to allow in the discovered offtargets. Note that Cas-OFFinder 3 is required to detect sites with bulges.",
        default=0,
        type=int,
    )
    makeguidefile_parser.add_argument(
        "-v",
        "--verbosity",
        help="Verbosity level of output to the console (1-4) 4 is the most verbose",
        type=int,
        default=3,
    )
    makeguidefile_parser.add_argument(
        "--debug",
        help="Print debug information",
        action="store_true",
    )

    parser.add_argument(
        "-v",
        "--verbosity",
        help="Verbosity level of output to the console (1-4) 4 is the most verbose",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--debug",
        help="Print debug information",
        action="store_true",
    )

    args = parser.parse_args()
    CRISPRessoShared.set_console_log_level(logger, args.verbosity, args.debug)

    if args.subcommand == "MakeGuideFile":
        output_folder = args.output_folder
        if output_folder is None:
            output_folder = "CRISPRessoSea_MakeGuideFileOutput"
        if not output_folder.endswith("/"):
            output_folder += "/"

        genome_file = args.genome_file
        if not os.path.isfile(genome_file):
            if os.path.isfile(args.genome_file+".fa"):
                genome_file = args.genome_file+".fa"
            if os.path.isfile(args.genome_file+".fasta"):
                genome_file = args.genome_file+".fasta"
        if not os.path.exists(genome_file):
            raise Exception('Cannot find genome fasta file at ' + args.genome_file)

        os.makedirs(output_folder, exist_ok=True)
        make_target_info_file(args.guide_seq, args.guide_name, args.pam, genome_file, args.max_mismatches, args.max_rna_bulges, args.max_dna_bulges, output_folder, args.file_prefix)
    
    elif args.subcommand == "Replot":
        if not os.path.isfile(args.reordered_stats_file):
            raise Exception('Reordered statistics file is not found at "' + args.reordered_stats_file + '". Use an aggregated_stats_all.txt file from a previous run as a template.')

        if not os.path.isfile(args.reordered_sample_file):
            raise Exception("Sample file not found at " + args.sample_file + "\n" + \
                            "Please provide a sample file with headers ['Name','fastq_r1','fastq_r2']")
        
        plot_group_order = None
        if args.plot_group_order is not None:
            plot_group_order = args.plot_group_order.split(",")
            for group in plot_group_order:
                if group not in list(pd.read_csv(args.sample_file, sep="\t")["group"].unique()):
                    raise Exception("In parsing group order (" + args.plot_group_order + ") group " + group + " not found in sample file")

        dot_plot_ylims = [None, None]
        if args.dot_plot_ylims is not None:
            try:
                dot_plot_ylim_vals = args.dot_plot_ylims.split(",")
                if len(dot_plot_ylim_vals) != 2:
                    raise Exception(
                        "Dot plot ylims must be a comma-separated list of two values"
                    )
                if dot_plot_ylim_vals[0] != "None":
                    dot_plot_ylims[0] = float(dot_plot_ylim_vals[0])
                if dot_plot_ylim_vals[1] != "None":
                    dot_plot_ylims[1] = float(dot_plot_ylim_vals[1])
            except Exception as e:
                raise Exception("Could not parse dot plot ylims: " + str(e))
        sig_method_parameters = None
        if args.sig_method_parameters is not None:
            sig_method_parameters = parse_sig_method_parameters(args.sig_method_parameters)

        replot(
            args.reordered_stats_file,
            args.reordered_sample_file,
            args.output_folder,
            args.file_prefix,
            name_column=args.name_column,
            fig_height=args.fig_height,
            fig_width=args.fig_width,
            seq_plot_ratio=args.seq_plot_ratio,
            plot_group_order=plot_group_order,
            dot_plot_ylims=dot_plot_ylims,
            heatmap_max_value=args.heatmap_max_value,
            heatmap_min_value=args.heatmap_min_value,
            sig_method_parameters=sig_method_parameters,
        )

    elif args.subcommand == "Process":
        if not os.path.isfile(args.sample_file):
            raise Exception("Sample file not found at " + args.sample_file + "\n" + \
                            "Please provide a sample file with headers ['Name','fastq_r1','fastq_r2']")
        # make sure sample file is in proper format
        _ = parse_sample_file(args.sample_file)

        if not os.path.isfile(args.target_file):
            raise Exception("Target file not found at " + args.target_file)
        # make sure target file is in proper format
        _ = parse_target_info(args.target_file)

        if args.gene_annotations is not None and not os.path.isfile(
            args.gene_annotations
        ):
            raise Exception(
                "Gene annotations file not found at " + args.gene_annotations
            )

        if not os.path.isfile(args.genome_file):
            raise Exception('Cannot find genome fasta file at ' + args.genome_file)
        index_files_exist = False
        if os.path.isfile(args.genome_file + ".1.bt2"):
            index_files_exist = True
        if os.path.isfile(args.genome_file + ".1.bt2l"):
            index_files_exist = True
        short_genome_file_name = '.'.join(args.genome_file.split(".")[0:-1])
        if os.path.isfile(short_genome_file_name + ".1.bt2"):
            index_files_exist = True
        if os.path.isfile(short_genome_file_name + ".1.bt2l"):
            index_files_exist = True
        if not index_files_exist:
            raise Exception(
                'Cannot find Bowtie2 index files for genome fasta file at ' + args.genome_file
            )

        n_processes_for_sea = 1
        if args.n_processes == "max":
            n_processes_for_sea = CRISPRessoMultiProcessing.get_max_processes()
        else:
            try:
                n_processes_for_sea = int(args.n_processes)
            except:
                raise Exception('Number of processors must be an integer or "max"')
        
        output_folder = args.output_folder
        if output_folder is None:
            output_folder = "CRISPRessoSea_output_on_" + os.path.basename(
                args.sample_file
            )
        if not output_folder.endswith("/"):
            output_folder += "/"
        os.makedirs(output_folder, exist_ok=True)

        plot_group_order = None
        if args.plot_group_order is not None:
            plot_group_order = args.plot_group_order.split(",")
            for group in plot_group_order:
                if group not in list(pd.read_csv(args.sample_file, sep="\t")["group"].unique()):
                    raise Exception("In parsing group order (" + args.plot_group_order + ") group " + group + " not found in sample file")

        sig_method_parameters = None
        if args.sig_method_parameters is not None:
            sig_method_parameters = parse_sig_method_parameters(args.sig_method_parameters)

        process_pools(
            args=args,
            sample_file=args.sample_file,
            target_file=args.target_file,
            genome_file=args.genome_file,
            output_folder=output_folder,
            gene_annotations=args.gene_annotations,
            n_processors=n_processes_for_sea,
            crispresso_quantification_window_center=args.crispresso_quantification_window_center,
            crispresso_quantification_window_size=args.crispresso_quantification_window_size,
            crispresso_base_editor_output=args.crispresso_base_editor_output,
            crispresso_default_min_aln_score=args.crispresso_default_min_aln_score,
            crispresso_plot_window_size=args.crispresso_plot_window_size,
            allow_unplaced_chrs=args.allow_unplaced_chrs,
            plot_only_complete_targets=args.plot_only_complete_targets,
            min_amplicon_coverage=args.min_amplicon_coverage,
            sort_based_on_mismatch=args.sort_based_on_mismatch,
            allow_target_match_to_other_region_loc=args.allow_target_match_to_other_region_loc,
            top_percent_cutoff=args.top_percent_cutoff,
            min_amplicon_len=args.min_amplicon_len,
            fail_on_pooled_fail=args.fail_on_pooled_fail,
            plot_group_order=plot_group_order,
            sig_method_parameters=sig_method_parameters,
        )
    else:
        raise Exception(
            'Please run with one of the following commands: \n' + 
            '"MakeGuideFile" (to make a target file from a given guide sequence)\n' +
            '"Process" (to process initial analysis) \n' + 
            '"Replot" (to replot finished analysis)'
        )
