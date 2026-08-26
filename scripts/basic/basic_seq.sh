#! /bin/bash
#Use ncov-sim as the environment (needs python 3.12)
#2e-6 is approximate sites per day evolution rate for SARS-CoV-2
phastSim \
    --seed 17 \
    --ref data/reference/wuhan.fasta \
    --tree results/basic/basic.nwk \
    --eteFormat 1 \
    --outpath results/basic/ \
    --createFasta \
    --scale 0.000002

#Use PairSNP to just get the identical sequences
pairsnp -s -d 0 results/basic/sars-cov-2_simulation_output.fasta > results/basic/idseq.tsv