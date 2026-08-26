#!/usr/bin/env bash
# Use ncov-sim as the environment (needs python 3.12)
# 2e-6 is approximate sites per day evolution rate for SARS-CoV-2

# Remove any previous RR TSV and simulation files
echo 'Hello!'
rm -f results/shutdown_threebody/rr.tsv
rm -f results/shutdown_threebody/transmission_pairs.tsv
rm -f results/shutdown_threebody/tree.nwk
rm -f results/shutdown_threebody/idseq.tsv
rm -f results/shutdown_threebody/meta.tsv
rm -f results/shutdown_threebody/sars-cov-2_simulation_output.fasta
rm -f results/shutdown_threebody/sars-cov-2_simulation_output.txt

for((i = 1; i <= 50; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    Rscript scripts/shutdown_threebody/shutdown_threebody_sim.R --seed $i

    # Note for phastSim we keep seed constant across all simulations
    phastSim \
        --seed 17 \
        --ref data/reference/wuhan.fasta \
        --tree results/shutdown_threebody/tree.nwk \
        --eteFormat 1 \
        --outpath results/shutdown_threebody/ \
        --createFasta \
        --scale 0.00001

    # Use PairSNP to just get the identical sequences
    pairsnp -s -t 6 -d 0 results/shutdown_threebody/sars-cov-2_simulation_output.fasta > results/shutdown_threebody/idseq.tsv

    # Analyze parent-child transmission pair sequence identity
    Rscript scripts/common/transmission_pair_analysis.R \
        --seed $i \
        --meta results/shutdown_threebody/meta.tsv \
        --pairs results/shutdown_threebody/idseq.tsv \
        --output results/shutdown_threebody/transmission_pairs.tsv
    Rscript scripts/shutdown_threebody/shutdown_threebody_RR.R --seed $i

    # Make sure the old iteration data files are not used for the next iteration if there is an error
    rm -f results/shutdown_threebody/tree.nwk
    rm -f results/shutdown_threebody/idseq.tsv
    rm -f results/shutdown_threebody/meta.tsv

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript scripts/shutdown_threebody/shutdown_threebody_plot.R
