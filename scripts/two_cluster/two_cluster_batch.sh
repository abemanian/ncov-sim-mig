#!/usr/bin/env bash
# Use ncov-sim as the environment (needs python 3.12)
# Four-deme two-cluster sweep: A-B and X-Y joined by the single A-X edge

SCENARIO="${1:-homogeneous}"
N_SEEDS="${2:-25}"
SUBSAMPLE_N="${3:-50000}"
INTERVENTION="${4:-0}"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 [SCENARIO] [N_SEEDS] [SUBSAMPLE_N] [INTERVENTION]"
    echo ""
    echo "SCENARIO: Intervention scenario (default: homogeneous)"
    echo "    none          - No intervention (null model)"
    echo "    homogeneous   - Staged movement schedule applied equally to all edges"
    echo "    heterogeneous - Same schedule, but A-X restoration lags by one phase"
    echo "    all           - All six arms (3 scenarios x intervention off/on) in"
    echo "                    sequence, then two_cluster_combined_plot.R. INTERVENTION"
    echo "                    is ignored. Per-arm logs go to results/two_cluster/logs/<date>/"
    echo ""
    echo "N_SEEDS:     Number of seeds to run (default: 25)"
    echo "SUBSAMPLE_N: Tips fed to pairsnp per replicate (default: 50000)"
    echo "INTERVENTION: 1 to reduce transmission in cluster 1 from t = 175"
    echo "              (sim --vaccine); outputs are prefixed <SCENARIO>_intervention"
    echo "              (default: 0)"
    exit 0
fi

# Full sweep: run every arm by calling this script once per arm. Sequential on
# purpose -- the arms share intermediate files in results/two_cluster/.
if [ "$SCENARIO" == "all" ]; then
    LOG_DIR="results/two_cluster/logs/$(date +%Y%m%d)"
    mkdir -p "$LOG_DIR"
    DRIVER_LOG="$LOG_DIR/driver.log"
    for scen in none homogeneous heterogeneous; do
        for iv in 0 1; do
            tag=$scen; [ "$iv" == "1" ] && tag="${scen}_intervention"
            echo "=== START $tag $(date '+%F %T') ===" >> "$DRIVER_LOG"
            bash "$0" "$scen" "$N_SEEDS" "$SUBSAMPLE_N" "$iv" > "$LOG_DIR/$tag.log" 2>&1
            echo "=== END $tag rc=$? $(date '+%F %T') ===" >> "$DRIVER_LOG"
        done
    done
    eval "$(conda shell.bash hook)"
    conda activate ncov-sim
    Rscript scripts/two_cluster/two_cluster_combined_plot.R >> "$DRIVER_LOG" 2>&1
    echo "=== ALL DONE $(date '+%F %T') ===" >> "$DRIVER_LOG"
    exit 0
fi

eval "$(conda shell.bash hook)"
conda activate ncov-sim

set -u

# The intervention arm gets its own output prefix (ARM) so it does not overwrite
# the matching non-intervention sweep. SCENARIO is still what the sim validates.
if [ "$INTERVENTION" == "1" ]; then
    ARM="${SCENARIO}_intervention"
    SIM_FLAG="--vaccine"
    ARM_FLAG="--intervention"
else
    ARM="$SCENARIO"
    SIM_FLAG=""
    ARM_FLAG=""
fi

# Clear per-scenario accumulators. Files are scenario-tagged in their contents,
# so a fresh sweep of one scenario should not inherit another's rows.
for f in rr rr_cluster rr_phase rr_cluster_phase rr_sliding rr_cluster_sliding \
         schedules rt_region; do
    rm -f "results/two_cluster/${ARM}_${f}.tsv"
done
rm -f results/two_cluster/tree.nwk results/two_cluster/tree_sub.nwk \
      results/two_cluster/subsample_ids.txt results/two_cluster/sim_info.tsv \
          results/two_cluster/idseq.tsv \
      results/two_cluster/meta.tsv results/two_cluster/schedule.tsv \
      results/two_cluster/subsampled.fasta \
      results/two_cluster/sars-cov-2_simulation_output.fasta \
      results/two_cluster/sars-cov-2_simulation_output.txt

echo "Scenario: $SCENARIO | intervention: $INTERVENTION | arm: $ARM | seeds: $N_SEEDS | subsample: $SUBSAMPLE_N"

for((i = 1; i <= N_SEEDS; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    if [ $i -eq 1 ]; then
        Rscript scripts/two_cluster/two_cluster_sim.R --seed $i --scenario $SCENARIO $SIM_FLAG --skip_tree --draw
    else
        Rscript scripts/two_cluster/two_cluster_sim.R --seed $i --scenario $SCENARIO $SIM_FLAG --skip_tree
    fi

    # Subsample BEFORE sequence simulation. The sim skips nosoi's full tree
    # (--skip_tree), which takes ~an hour for a 250k-host outbreak; instead the
    # tips are picked from meta.tsv (time-stratified) and the tree of just those
    # tips is built straight from the transmission chain. phastSim then only
    # writes SUBSAMPLE_N sequences. Equivalence with pruning the full nosoi tree
    # is checked by scripts/common/validate_sampled_tree.R.
    python scripts/common/subsample_fasta.py \
        --ids-out results/two_cluster/subsample_ids.txt \
        --meta results/two_cluster/meta.tsv \
        --n $SUBSAMPLE_N \
        --seed 17 \
        --bin-width 10

    Rscript scripts/common/build_sampled_tree.R \
        --meta results/two_cluster/meta.tsv \
        --ids results/two_cluster/subsample_ids.txt \
        --total-time "$(tail -n 1 results/two_cluster/sim_info.tsv | cut -f 2)" \
        --output results/two_cluster/tree_sub.nwk

    # phastSim seed is held constant so stochasticity stays in the epi model
    phastSim \
        --seed 17 \
        --ref data/reference/wuhan.fasta \
        --tree results/two_cluster/tree_sub.nwk \
        --eteFormat 1 \
        --outpath results/two_cluster/ \
        --createFasta \
        --scale 0.000002

    mv results/two_cluster/sars-cov-2_simulation_output.fasta \
       results/two_cluster/subsampled.fasta

    pairsnp -s -t 8 -d 0 results/two_cluster/subsampled.fasta > results/two_cluster/idseq.tsv

    Rscript scripts/two_cluster/two_cluster_RR.R --seed $i --scenario $SCENARIO $ARM_FLAG
    Rscript scripts/two_cluster/two_cluster_reproduction_region.R --seed $i --scenario $SCENARIO $ARM_FLAG

    rm -f results/two_cluster/tree.nwk results/two_cluster/tree_sub.nwk \
          results/two_cluster/subsample_ids.txt results/two_cluster/sim_info.tsv \
          results/two_cluster/idseq.tsv \
          results/two_cluster/meta.tsv results/two_cluster/schedule.tsv \
          results/two_cluster/subsampled.fasta \
          results/two_cluster/sars-cov-2_simulation_output.fasta

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript scripts/two_cluster/two_cluster_plot.R --scenario $SCENARIO $ARM_FLAG
