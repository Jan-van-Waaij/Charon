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

#=
numberofindividualsinthesedimentvec = 6:6
coverages = [50, 40, 39:-1:20..., 15, 10]
allcombinations = [(numind, cov) for numind in numberofindividualsinthesedimentvec for cov in coverages]
=#

numberofindividualsinthesedimentvec = 1:8
coverages = 100:-1:1

@assert issorted(coverages, rev=true)
@assert allunique(coverages)

for numberofindividualsinthesediment in numberofindividualsinthesedimentvec
    for index in eachindex(coverages)[1:end-1]
        oldcoverage = coverages[index]
        newcoverage = coverages[index+1]
        @show oldcoverage newcoverage
        for index in eachindex(numberofindividualsinthesedimentvec)
            readfileoldcoverage = string("data/reads/subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", oldcoverage, "_reads_with_errors.txt.gz")
            readfilenewcoverage = string("data/reads/subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", newcoverage, "_reads_with_errors.txt.gz")
            comments = """
            subsample from $readfileoldcoverage
            new coverage: $newcoverage
            """
            subsample(readfileoldcoverage, readfilenewcoverage, oldcoverage, newcoverage, comments)
        end
    end 
end 


#= tests 
for numberofindividualsinthesediment in numberofindividualsinthesedimentvec
    for index in eachindex(coverages)[1:end-1]
        oldcoverage = coverages[index]
        newcoverage = coverages[index+1]
        for index in eachindex(numberofindividualsinthesedimentvec)
            readfileoldcoverage = string("data/reads/subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", oldcoverage, "_reads_with_errors.txt.gz")
            readfilenewcoverage = string("data/reads/subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", newcoverage, "_reads_with_errors.txt.gz")
            oldder, oldcov = parsereadfile(readfileoldcoverage)
            newder, newcov = parsereadfile(readfilenewcoverage)
            @test newder ≤ oldder 
            @test length(unique(oldcov)) == 1 
            @test length(unique(newcov)) == 1 
            @test oldcov[1] == oldcoverage
            @test newcov[1] == newcoverage

        end
    end 
end 
=#
