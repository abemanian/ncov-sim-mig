#! /bin/bash
#Use ncov-sim as the environment (needs python 3.12)
#2e-6 is approximate sites per day evolution rate for SARS-CoV-2
phastSim \
    --seed 17 \
    --ref wuhan.fasta \
    --tree shutdown_three/tree.nwk \
    --eteFormat 1 \
    --outpath shutdown_three/ \
    --createFasta \
    --scale 0.00001

#Use PairSNP to just get the identical sequences
pairsnp -s -d 0 shutdown_three/sars-cov-2_simulation_output.fasta > shutdown_three/idseq.tsv