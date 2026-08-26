#! /bin/bash
# Use ncov-sim as the environment (needs python 3.12)
# Single-replicate sequence evolution for the two-cluster simulation
SUBSAMPLE_N="${1:-8000}"

phastSim \
    --seed 17 \
    --ref data/reference/wuhan.fasta \
    --tree results/two_cluster/tree.nwk \
    --eteFormat 1 \
    --outpath results/two_cluster/ \
    --createFasta \
    --scale 0.000002

# pairsnp is O(n^2); subsample tips to keep the pairwise step tractable
python scripts/common/subsample_fasta.py \
    --input results/two_cluster/sars-cov-2_simulation_output.fasta \
    --output results/two_cluster/subsampled.fasta \
    --n $SUBSAMPLE_N \
    --seed 17 \
    --meta results/two_cluster/meta.tsv \
    --bin-width 10

pairsnp -s -d 0 results/two_cluster/subsampled.fasta > results/two_cluster/idseq.tsv
