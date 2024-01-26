#!/bin/bash

# Run msprime_sim_v1.py and pipe the output to bgzip and save as sim_v1_vcf.gz
nice -19 python3 msprime_isolates.py ../demes/OOA_KAR_extension_splits.yaml mutation_ancestry_positions.tsv | bgzip -c > isoaltes_4pop_1kindi.vcf.gz
