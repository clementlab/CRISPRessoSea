# set --top_percent_cutoff to 80% because we have a lot of regions and they wouldn't be picked up with the default 20% cutoff
CRISPRessoSea Process -s samples.small.txt -t CRISPRessoSea_MakeGuideFileOutput/CRISPRessoSea.guide_info.txt --genome_file test_genome.fa --plot_only_complete_targets --top_percent_cutoff 0.8 --n_processes 24 --verbosity 4
