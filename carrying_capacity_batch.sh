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
rm -f carrying_three/rr.tsv
rm -f carrying_three/rr_sliding.tsv
rm -f carrying_three/transmission_pairs.tsv
rm -f carrying_three/tree.nwk
rm -f carrying_three/idseq.tsv
rm -f carrying_three/meta.tsv
rm -f carrying_three/sars-cov-2_simulation_output.fasta
rm -f carrying_three/sars-cov-2_simulation_output.txt

for((i = 1; i <= 50; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    #Draw the initial simulation diagnostic plots for the first plot
    if [ $i -eq 1 ]; then
        Rscript carrying_capacity_sim.R --seed $i --scenario $SCENARIO --draw
    else
        Rscript carrying_capacity_sim.R --seed $i --scenario $SCENARIO
    fi

    # Note for phastSim we keep seed constant across all simulations
    phastSim \
        --seed 17 \
        --ref wuhan.fasta \
        --tree carrying_three/tree.nwk \
        --eteFormat 1 \
        --outpath carrying_three/ \
        --createFasta \
        --scale 0.000002

    # Use PairSNP to just get the identical sequences
    pairsnp -s -t 16 -d 0 carrying_three/sars-cov-2_simulation_output.fasta > carrying_three/idseq.tsv

    # Analyze parent-child transmission pair sequence identity
    Rscript transmission_pair_analysis.R \
        --seed $i \
        --meta carrying_three/meta.tsv \
        --pairs carrying_three/idseq.tsv \
        --output carrying_three/transmission_pairs.tsv

    Rscript carrying_capacity_RR.R --seed $i

    # Make sure the old iteration data files are not used for the next iteration if there is an error
    rm -f carrying_three/tree.nwk
    rm -f carrying_three/idseq.tsv
    rm -f carrying_three/meta.tsv

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript carrying_capacity_plot.R --scenario $SCENARIO