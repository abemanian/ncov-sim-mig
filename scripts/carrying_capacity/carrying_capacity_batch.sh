#!/usr/bin/env bash
# Use ncov-sim as the environment (needs python 3.12)
# 2e-6 is approximate sites per day evolution rate for SARS-CoV-2

# Parse command-line arguments
SCENARIO="${1:-basic}"

# Show usage if help is requested
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 [SCENARIO]"
    echo ""
    echo "SCENARIO: Simulation scenario to run (default: basic)"
    echo "  Options:"
    echo "    basic              - Basic transmission model"
    echo "    seasonal_cyclic    - Seasonal transmission with cyclical force of infection"
    echo "    seasonal_threshold - Seasonal transmission where force of infection is decreased for a specific time window"
    echo "    isolated           - Model with isolation interventions (deme A locked down at given time point)"
    echo "    immune             - Model with immunization interventions (80% of deme A immunized at given time point"
    exit 0
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ncov-sim

# Remove any previous RR TSV and simulation files
rm -f results/carrying_capacity/rr.tsv
rm -f results/carrying_capacity/rr_sliding.tsv
rm -f results/carrying_capacity/transmission_pairs.tsv
rm -f results/carrying_capacity/tree.nwk
rm -f results/carrying_capacity/idseq.tsv
rm -f results/carrying_capacity/meta.tsv
rm -f results/carrying_capacity/sars-cov-2_simulation_output.fasta
rm -f results/carrying_capacity/sars-cov-2_simulation_output.txt

for((i = 1; i <= 50; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    #Draw the initial simulation diagnostic plots for the first plot
    if [ $i -eq 1 ]; then
        Rscript scripts/carrying_capacity/carrying_capacity_sim.R --seed $i --scenario $SCENARIO --draw
    else
        Rscript scripts/carrying_capacity/carrying_capacity_sim.R --seed $i --scenario $SCENARIO
    fi

    # Note for phastSim we keep seed constant across all simulations
    phastSim \
        --seed 17 \
        --ref data/reference/wuhan.fasta \
        --tree results/carrying_capacity/tree.nwk \
        --eteFormat 1 \
        --outpath results/carrying_capacity/ \
        --createFasta \
        --scale 0.000002

    # Use PairSNP to just get the identical sequences
    pairsnp -s -t 16 -d 0 results/carrying_capacity/sars-cov-2_simulation_output.fasta > results/carrying_capacity/idseq.tsv

    # Analyze parent-child transmission pair sequence identity
    Rscript scripts/common/transmission_pair_analysis.R \
        --seed $i \
        --meta results/carrying_capacity/meta.tsv \
        --pairs results/carrying_capacity/idseq.tsv \
        --output results/carrying_capacity/transmission_pairs.tsv

    Rscript scripts/carrying_capacity/carrying_capacity_RR.R --seed $i

    # Make sure the old iteration data files are not used for the next iteration if there is an error
    rm -f results/carrying_capacity/tree.nwk
    rm -f results/carrying_capacity/idseq.tsv
    rm -f results/carrying_capacity/meta.tsv

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript scripts/carrying_capacity/carrying_capacity_plot.R --scenario $SCENARIO