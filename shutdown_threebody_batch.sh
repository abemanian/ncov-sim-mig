#!/usr/bin/env bash
# Use ncov-sim as the environment (needs python 3.12)
# 2e-6 is approximate sites per day evolution rate for SARS-CoV-2

# Remove any previous RR TSV and simulation files
rm -f shutdown_three/rr.tsv
rm -f shutdown_three/tree.nwk
rm -f shutdown_three/idseq.tsv
rm -f shutdown_three/meta.tsv

for((i = 1; i <= 50; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    Rscript shutdown_threebody_sim.R --seed $i

    # Note for phastSim we keep seed constant across all simulations
    phastSim \
        --seed 17 \
        --ref wuhan.fasta \
        --tree shutdown_three/tree.nwk \
        --eteFormat 1 \
        --outpath shutdown_three/ \
        --createFasta \
        --scale 0.00001

    # Use PairSNP to just get the identical sequences
    pairsnp -s -t 6 -d 0 shutdown_three/sars-cov-2_simulation_output.fasta > shutdown_three/idseq.tsv

    Rscript shutdown_threebody_RR.R --seed $i

    # Make sure the old iteration data files are not used for the next iteration if there is an error
    rm -f shutdown_three/tree.nwk
    rm -f shutdown_three/idseq.tsv
    rm -f shutdown_three/meta.tsv

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript shutdown_threebody_plot.R
