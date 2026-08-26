#!/usr/bin/env bash
# Use ncov-sim as the environment (needs python 3.12)
# Six-deme two-cluster sweep: A-B-C and X-Y-Z joined by the single A-X edge

SCENARIO="${1:-homogeneous}"
N_SEEDS="${2:-25}"
SUBSAMPLE_N="${3:-8000}"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 [SCENARIO] [N_SEEDS] [SUBSAMPLE_N]"
    echo ""
    echo "SCENARIO: Intervention scenario (default: homogeneous)"
    echo "    none          - No intervention (null model)"
    echo "    homogeneous   - Staged movement schedule applied equally to all edges"
    echo "    heterogeneous - Same schedule, but A-X restoration lags by one phase"
    echo ""
    echo "N_SEEDS:     Number of seeds to run (default: 25)"
    echo "SUBSAMPLE_N: Tips fed to pairsnp per replicate (default: 8000)"
    exit 0
fi

eval "$(conda shell.bash hook)"
conda activate ncov-sim

set -u

# Clear per-scenario accumulators. Files are scenario-tagged in their contents,
# so a fresh sweep of one scenario should not inherit another's rows.
for f in rr rr_cluster rr_phase rr_cluster_phase rr_sliding rr_cluster_sliding \
         schedules; do
    rm -f "results/two_cluster/${SCENARIO}_${f}.tsv"
done
rm -f "results/two_cluster/${SCENARIO}_transmission_pairs.tsv"
rm -f results/two_cluster/tree.nwk results/two_cluster/idseq.tsv \
      results/two_cluster/meta.tsv results/two_cluster/schedule.tsv \
      results/two_cluster/subsampled.fasta \
      results/two_cluster/sars-cov-2_simulation_output.fasta \
      results/two_cluster/sars-cov-2_simulation_output.txt

echo "Scenario: $SCENARIO | seeds: $N_SEEDS | subsample: $SUBSAMPLE_N"

for((i = 1; i <= N_SEEDS; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    if [ $i -eq 1 ]; then
        Rscript scripts/two_cluster/two_cluster_sim.R --seed $i --scenario $SCENARIO --draw
    else
        Rscript scripts/two_cluster/two_cluster_sim.R --seed $i --scenario $SCENARIO
    fi

    # phastSim seed is held constant so stochasticity stays in the epi model
    phastSim \
        --seed 17 \
        --ref data/reference/wuhan.fasta \
        --tree results/two_cluster/tree.nwk \
        --eteFormat 1 \
        --outpath results/two_cluster/ \
        --createFasta \
        --scale 0.000002

    python scripts/common/subsample_fasta.py \
        --input results/two_cluster/sars-cov-2_simulation_output.fasta \
        --output results/two_cluster/subsampled.fasta \
        --n $SUBSAMPLE_N \
        --seed 17 \
        --meta results/two_cluster/meta.tsv \
        --bin-width 10

    pairsnp -s -t 8 -d 0 results/two_cluster/subsampled.fasta > results/two_cluster/idseq.tsv

    Rscript scripts/common/transmission_pair_analysis.R \
        --seed $i \
        --meta results/two_cluster/meta.tsv \
        --pairs results/two_cluster/idseq.tsv \
        --output results/two_cluster/${SCENARIO}_transmission_pairs.tsv

    Rscript scripts/two_cluster/two_cluster_RR.R --seed $i --scenario $SCENARIO

    rm -f results/two_cluster/tree.nwk results/two_cluster/idseq.tsv \
          results/two_cluster/meta.tsv results/two_cluster/schedule.tsv \
          results/two_cluster/subsampled.fasta \
          results/two_cluster/sars-cov-2_simulation_output.fasta

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript scripts/two_cluster/two_cluster_plot.R --scenario $SCENARIO
