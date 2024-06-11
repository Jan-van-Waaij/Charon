# simreadssubsampling.jl
# Author: Jan van Waaij
# Created: 10 June 2024
# Last update: 10 June 2024
#
# Purpose:
# Simulate reads for the numberofindividualsinthesediment individuals in the sediment, 
# who are randomly chosen from the ancient population in filename.vcf.gz. The coverage is set to 100. 
#
# Installation:
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open file as julia path/to/simreadssubsampling.jl
#
# Output:
# output will be saved in data/reads/subsampling_number_of_individuals=$numberofindividualsinthesediment_coverage=100_reads_with_errors.txt.gz.
# The output will have a few lines of comments, each of which preceded by #, and then each line has two integers, of which the first is the number
# of derived reads and the second is the total number of reads (a.k.a the coverage).

using Data, GZip, StatsBase

infoblock = 9
ancientpopulation = infoblock .+ (1:5000)

vcffile = "data/vcffiles_demes/two_demes_with_demes_tauC=0.25tauA=0.25numgen=1000effpop=2000_12.vcf.gz"

numberofindividualsinthesedimentvec = 1:8 # [1, 2, 4, 8]
ϵ = 0.003

randomindividuals = similar(numberofindividualsinthesedimentvec, Vector{Int})

@assert issorted(numberofindividualsinthesedimentvec)
randomindividuals[end] = sample(ancientpopulation, numberofindividualsinthesedimentvec[end], replace=false, ordered = true)
for k in eachindex(numberofindividualsinthesedimentvec, randomindividuals)[end-1:-1:begin]
    randomindividuals[k] = sample(randomindividuals[k+1], numberofindividualsinthesedimentvec[k], replace=false, ordered=true)
end 

@assert length(randomindividuals) == length(numberofindividualsinthesedimentvec)
for k in eachindex(numberofindividualsinthesedimentvec)
    @assert allunique(randomindividuals[k])
    @assert length(randomindividuals[k]) == numberofindividualsinthesedimentvec[k] 
    if k != lastindex(randomindividuals)
        @assert randomindividuals[k] ⊆ randomindividuals[k+1]
    end 
end

#=
numberofindividualsinthesedimentvec = [2, 4, 8]
for _ in 1:100
    randomindividuals = similar(numberofindividualsinthesedimentvec, Vector{Int})
    @assert issorted(numberofindividualsinthesedimentvec)
    randomindividuals[end] = sample(ancientpopulation, numberofindividualsinthesedimentvec[end], replace=false, ordered = true)
    for k in eachindex(numberofindividualsinthesedimentvec)[end-1:-1:begin]
        randomindividuals[k] = sample(randomindividuals[k+1], numberofindividualsinthesedimentvec[k], replace=false, ordered=true)
    end 
    @test length(randomindividuals) == length(numberofindividualsinthesedimentvec)
    for k in eachindex(numberofindividualsinthesedimentvec)
        @test allunique(randomindividuals[k])
        @test length(randomindividuals[k]) == numberofindividualsinthesedimentvec[k] 
        if k != lastindex(numberofindividualsinthesedimentvec)
            @test randomindividuals[k] ⊆ randomindividuals[k+1]
        end 
    end 
end
=#

for index in eachindex(numberofindividualsinthesedimentvec, randomindividuals)
    # I choose numberofindividualsinthesediment individuals from the ancient population.
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
