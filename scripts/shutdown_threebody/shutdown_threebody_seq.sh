#! /bin/bash
#Use ncov-sim as the environment (needs python 3.12)
#2e-6 is approximate sites per day evolution rate for SARS-CoV-2
phastSim \
    --seed 17 \
    --ref data/reference/wuhan.fasta \
    --tree results/shutdown_threebody/tree.nwk \
    --eteFormat 1 \
    --outpath results/shutdown_threebody/ \
    --createFasta \
    --scale 0.00001

#Use PairSNP to just get the identical sequences
pairsnp -s -d 0 results/shutdown_threebody/sars-cov-2_simulation_output.fasta > results/shutdown_threebody/idseq.tsv