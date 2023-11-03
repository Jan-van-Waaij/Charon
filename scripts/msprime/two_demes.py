import sys
import msprime
import numpy as np
import scipy.stats
import matplotlib.pyplot as pyplot
import math
import random
import subprocess

#POP1,POP2 = 0,1 
#split_time=14;


def two_demes( effectivePop1,effectivePop2, joinTimeYears, randomseed ):
    #infile = "hapmap/genetic_map_GRCh37_chr1.txt"
    # Read in the recombination map using the read_hapmap method,
    #recomb_map = msprime.RecombinationMap.read_hapmap(infile)

    # Times are provided in years, so we convert into generations.
    generation_time = 25
    #joinTime        = 21200 / generation_time
    joinTime        = joinTimeYears / generation_time


    #we say that there are 2 populations with 10k individuals as sample size (not the Ne), it has 5k individuals in "diploid"
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=10000, initial_size=effectivePop1),#target population            population 0 in msprime
        msprime.PopulationConfiguration(sample_size=10000, initial_size=effectivePop2) #allele frequency population  population 1 in msprime
    ]

    #no one goes anywhere before jointime
    migration_matrix = [
        [0,0],
        [0,0],
    ]

    #
    demographic_events = [
        # CEU and CHB merge into B with rate changes at joinTime
        msprime.MassMigration(              time=joinTime, source=1, destination=0, proportion=1.0), #everyone in population 1 is now in population 0 (join) 
        msprime.MigrationRateChange(        time=joinTime, rate=0),
        msprime.MigrationRateChange(        time=joinTime, rate=0, matrix_index=(0, 1)),  
        msprime.MigrationRateChange(        time=joinTime, rate=0, matrix_index=(1, 0)),
        #                                                    once it is joined, the Ne is the max of both
        msprime.PopulationParametersChange( time=joinTime, initial_size=max(effectivePop1,effectivePop2), growth_rate=0, population_id=1), #it is now the max of both 
    ]      

    # Now we run the simulation for this map. We assume Ne=10^4
    # and have a sample of 100 individuals
    tree_sequence = msprime.simulate(
        #sample_size=16,
        length = 10e6,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        Ne=max(effectivePop1,effectivePop2),
        mutation_rate=2e-8,
        random_seed=10,
        recombination_rate=1e-8
        #recombination_map=recomb_map
    )
    
    tree_sequence.dump( "two_demes_"+str(effectivePop1)+"_"+str(effectivePop2)+"_"+str(joinTimeYears)+"_"+str(randomseed)+".tree");# ,  True);


    vcffilename="two_demes_"+str(effectivePop1)+"_"+str(effectivePop2)+"_"+str(joinTimeYears)+"_"+str(randomseed)+".vcf";
    with open(vcffilename, "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, 2)
    subprocess.call("/home/ctools/bin/bgzip -f "+vcffilename, shell=True)

#           effectivePop1(target),effectivePop2(allelefreq), joinTimeYears,    seed
two_demes( int(sys.argv[1]),      int(sys.argv[2]),          int(sys.argv[3]), int(sys.argv[4]) );

