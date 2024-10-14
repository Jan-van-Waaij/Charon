# simreadsMadsmodel.jl
# Author: Jan van Waaij
# Created: 26 February 2024
# Last update: 26 February 2024
#
# Purpose:
# Simulate reads from Mads model 
#
# Installation:
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open file as julia path/to/simreadsMadsmodel.jl
#
# Output:
# output will be saved in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads/MadsModel_population_coverage_reads.txt.gz
# The output will have a few lines of comments, each of which preceded by #, and then each line has two integers, of which the first is the number
# of derived reads and the second is the total number of reads (a.k.a the coverage).

using Data, GZip, StatsBase, Base.Threads

infoblock = 9
ancientpopulation = infoblock .+ (1:5000)

vcffile = "data/vcffilesMads/isoaltes_4pop_1kindi.vcf.gz"
readfilefirstpart = "data/reads/Mads_model_pop_"

populations = Dict(:CEU => 9 .+ (1:2000), :YRI => 9 .+ (2001:4000), :CHB => 9 .+ (4001:6000), :KAR => 9 .+ (6001:8000), :CEUS => 9 .+ (8001:10_000), :YRIS => 9 .+ (10_001:12_000), :CHBS => 9 .+ (12_001:14_000), :KARS => 9 .+ (14_001:16_000))

for pop in values(populations)
    @assert length(pop) == 2000
end 

numberofindividualsinthesedimentvec = 1:8
ϵ = 0.003
coverage = 100

populationnames = [:CHBS, :YRIS, :CEUS, :KARS]
allcombinations = [(population, numberofindividualsinthesediment) for population in populationnames for numberofindividualsinthesediment in numberofindividualsinthesedimentvec] # 32 elements. 
@threads for (population, numberofindividualsinthesediment) in allcombinations
    # I choose numberofindividualsinthesediment individuals from the ancient population. 
    individualsinthesediment = sample(populations[population], numberofindividualsinthesediment, 
                                        replace=false, ordered=true)

    # I open the file where I while store the reads. 
    local openedvcffile = GZip.open(vcffile)
    # skip the header. 
    skiplines(openedvcffile, 6)
    readfile = string(readfilefirstpart, population, "_num_ind_",
        numberofindividualsinthesediment, "_cov_", coverage, "_reads_with_errors.txt.gz")
    openedreadfile = GZip.open(readfile, "w")

    comments = """
    Each line contains the number of derived reads and the total number of reads, seperated by a tab.  
    original vcf file: $vcffile
    Individuals in the sediment: $individualsinthesediment
    Total number of individuals in the sediment: $numberofindividualsinthesediment
    coverage: $coverage
    error rate: $ϵ"""

    reads(openedvcffile, openedreadfile, individualsinthesediment, coverage, comments, ϵ)
    
    close(openedreadfile)
    close(openedvcffile)  
end 