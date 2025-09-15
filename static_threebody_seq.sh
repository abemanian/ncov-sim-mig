#! /bin/bash
#Use ncov-sim as the environment (needs python 3.12)
#2e-6 is approximate sites per day evolution rate for SARS-CoV-2
phastSim \
    --seed 17 \
    --ref wuhan.fasta \
    --tree static_three/static_threebody.nwk \
    --eteFormat 1 \
    --outpath static_three/ \
    --createFasta \
    --scale 0.000002

#Use PairSNP to just get the identical sequences
pairsnp -s -d 0 static_three/sars-cov-2_simulation_output.fasta > static_three/idseq.tsv