#!/usr/bin/env bash
# Use ncov-sim as the environment (needs python 3.12)
# 2e-6 is approximate sites per day evolution rate for SARS-CoV-2

# Remove any previous RR TSV and simulation files
rm -f carrying_three/rr.tsv
rm -f carrying_three/rr_sliding.tsv
rm -f carrying_three/tree.nwk
rm -f carrying_three/idseq.tsv
rm -f carrying_three/meta.tsv
rm -f carrying_three/sars-cov-2_simulation_output.fasta
rm -f carrying_three/sars-cov-2_simulation_output.txt

for((i = 1; i <= 50; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    Rscript carrying_capacity_sim.R --seed $i

    # Note for phastSim we keep seed constant across all simulations
    phastSim \
        --seed 17 \
        --ref wuhan.fasta \
        --tree carrying_three/tree.nwk \
        --eteFormat 1 \
        --outpath carrying_three/ \
        --createFasta \
        --scale 0.00001

    # Use PairSNP to just get the identical sequences
    pairsnp -s -t 14 -d 0 carrying_three/sars-cov-2_simulation_output.fasta > carrying_three/idseq.tsv

    Rscript carrying_capacity_RR.R --seed $i

    # Make sure the old iteration data files are not used for the next iteration if there is an error
    rm -f carrying_three/tree.nwk
    rm -f carrying_three/idseq.tsv
    rm -f carrying_three/meta.tsv

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript carrying_capacity_plot.R