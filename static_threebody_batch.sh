#!/usr/bin/env bash
#Use ncov-sim as the environment (needs python 3.12)
#2e-6 is approximate sites per day evolution rate for SARS-CoV-2

#Remove any previous RR TSV and simulation files
rm -f static_three/rr.tsv
rm -f static_three/tree.nwk
rm -f static_three/idseq.tsv
rm -f static_three/meta.tsv

for((i = 1; i <= 50; i ++)); do
    echo "Seed: $i"
    SECONDS=0

    Rscript static_threebody_sim.R --seed $i

    #Note for phastSim we keep seed constant across all simulations
    phastSim \
        --seed 17 \
        --ref wuhan.fasta \
        --tree static_three/tree.nwk \
        --eteFormat 1 \
        --outpath static_three/ \
        --createFasta \
        --scale 0.00001

    #Use PairSNP to just get the identical sequences
    pairsnp -s -t 6 -d 0 static_three/sars-cov-2_simulation_output.fasta > static_three/idseq.tsv

    Rscript static_threebody_RR.R --seed $i

    #Make sure the old iteration data files are not used for the next iteration if there is an error
    rm -f static_three/tree.nwk
    rm -f static_three/idseq.tsv
    rm -f static_three/meta.tsv

    echo "Finished Seed $i in ${SECONDS} seconds"
done

Rscript static_threebody_plot.R