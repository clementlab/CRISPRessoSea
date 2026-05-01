# Use 50 processes for the CRISPRessoSea run.
# Set top_percent_cutoff to 1 because these are simulated reads,
# so we want to include all aligned sites rather than only the top 20% of aligned locations (which would be expected if this were real data with noisy PCR amplification).
CRISPRessoSea Process -o 02_run_output \
    -t 01_make_simulations_output/guides.txt \
    -s 01_make_simulations_output/samples.txt \
    -x 01_make_simulations_output/reference.fasta \
    -p 50 \
    --top_percent_cutoff 1
