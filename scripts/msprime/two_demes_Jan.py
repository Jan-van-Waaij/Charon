import sys # for command line arguments. 
import msprime 
#import numpy as np
#import scipy.stats
#import matplotlib.pyplot as pyplot
#import math
#import random
import subprocess # to call gzip from within Python. 

def two_demes(NeAncient, NeAnchor, split_time_in_generations, random_seed=None):
    demography = msprime.Demography()
    demography.add_population(name="ancient", initial_size=NeAncient, growth_rate=0.0)
    demography.add_population(name="anchor", initial_size=NeAnchor, growth_rate=0.0)
    demography.add_population(name="common_ancestor", initial_size=max(NeAncient, NeAnchor), growth_rate=0.0)
    demography.add_population_split(time=split_time_in_generations, derived=["ancient", "anchor"], ancestral="common_ancestor")
    ts = msprime.sim_ancestry(samples={"ancient":5_000, "anchor":5_000}, demography=demography, random_seed=random_seed, recombination_rate=1e-8, sequence_length=10e6)
    mts = msprime.sim_mutations(ts, rate=2e-8)
    mts.dump( "two_demes_time_in_generations_"+str(NeAncient)+"_"+str(NeAnchor)+"_"+str(split_time_in_generations)+"_"+str(random_seed)+".tree");# ,  True);
    vcffilename="two_demes_time_in_generations_"+str(NeAncient)+"_"+str(NeAnchor)+"_"+str(split_time_in_generations)+"_"+str(random_seed)+".vcf";
    with open(vcffilename, "w") as vcf_file:
        mts.write_vcf(vcf_file)
    subprocess.call("/home/ctools/bin/bgzip -f "+vcffilename, shell=True)

#     NeAncient, NeAnchor, split_time_in_generations, random_seed=None
two_demes( int(sys.argv[1]),      int(sys.argv[2]),          int(sys.argv[3]), int(sys.argv[4]) );

# for tauC = tauA = 2, take NeAncient=NeAnchor=1000 and number_generations=4000

# py scripts/msprime/two_demes_Jan.py 1000 1000 4000 1234
# two_demes( int(sys.argv[1]),      int(sys.argv[2]),          int(sys.argv[3]), int(sys.argv[4]) );