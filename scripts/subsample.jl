# subsample.jl
# Author: Jan van Waaij
# Created: 10 June 2024
# Last update: 10 June 2024
#
# Purpose:
# subsample the reads to lower coverage.
#
# Installation:
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open file as julia path/to/subsample.jl
#
# Output:
# output will be saved in data/reads/subsampling_number_of_individuals=$numberofindividualsinthesediment_coverage=$coverage_reads_with_errors.txt.gz.
# The output will have a few lines of comments, each of which preceded by #, and then each line has two integers, of which the first is the number
# of derived reads and the second is the total number of reads (a.k.a the coverage).

using Data, GZip, StatsBase

numberofindividualsinthesedimentvec = 1:8 # [1, 2, 4, 8]

for index in eachindex(numberofindividualsinthesedimentvec)
    numberofindividualsinthesediment = numberofindividualsinthesedimentvec[index]
    individualsinthesediment = randomindividuals[index] 
    @assert length(randomindividuals[index]) == numberofindividualsinthesedimentvec[index]

    # I open the file where I while store the reads. 
    local openedvcffile = GZip.open(vcffile)
    # skip the header. 
    skiplines(openedvcffile, 6)
    readfile = string("data/reads/subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=100_reads_with_errors.txt.gz")
    openedreadfile = GZip.open(readfile, "w")

    comments = """
    Each line contains the number of derived reads and the total number of reads, seperated by a tab.  
    original vcf file: $vcffile
    Individuals in the sediment: $individualsinthesediment
    Total number of individuals in the sediment: $numberofindividualsinthesediment
    coverage: 100
    error rate: $ϵ"""

    reads(openedvcffile, openedreadfile, individualsinthesediment, 100, comments, ϵ)
    close(openedreadfile)
    close(openedvcffile)
end  
