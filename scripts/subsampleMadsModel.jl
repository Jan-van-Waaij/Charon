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

vcffile = "data/vcffilesMads/isoaltes_4pop_1kindi.vcf.gz"
readfilefirstpart = "data/reads/Mads_model_pop_"

populationnames = [:CHBS, :YRIS, :CEUS, :KARS]

for population in populationnames
    for numberofindividualsinthesediment in 1:8
        for oldcoverage in 100:-1:2
            newcoverage = oldcoverage - 1 
            readfileoldcoverage = string(readfilefirstpart, population, "_num_ind_",
            numberofindividualsinthesediment, "_cov_", oldcoverage, "_reads_with_errors.txt.gz")
            readfilenewcoverage = string(readfilefirstpart, population, "_num_ind_",
            numberofindividualsinthesediment, "_cov_", newcoverage, "_reads_with_errors.txt.gz")
            comments = """
            Each line contains the number of derived reads and the total number of reads, seperated by a tab.  
            original vcf file: $vcffile
            Total number of individuals in the sediment: $numberofindividualsinthesediment
            coverage: $newcoverage"""
            subsample(readfileoldcoverage, readfilenewcoverage, oldcoverage, newcoverage, comments)
        end 
    end 
end 
