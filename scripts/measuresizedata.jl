# runmcmconallfiles.jl
# Author: Jan van Waaij
# Created: 26 February 2024
# Last update: 27 February 2024
#
# Purpose:
# Run MCMC with the files simulated by Mads.
#
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia -O3 --threads=auto path/to/runtraditionalmcmcMadsmodel.jl
#

# Als compute01 ongebruikt is, run het daar, het heeft meer dan 128 CPUs. 


using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      StatsBase 

# /home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind16_basecount_IBS_ind16.gz
# let op de coverage is niet constant. 
anchorpops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
anchorpops = sample(anchorpops, 3, replace=false)
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, truenumind) for anchorpop in anchorpops for truenumind in truenumindvec] # 104 elementen. 4*104=416 threads



results = Tuple{Vector{Int64}, Vector{Int64}, Vector{Float64}}[]
for (anchorpop, truenumind) in allcombinations  
    ancientpop = "IBS"  
    frequencyfile = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind$(truenumind)_freq_$(anchorpop).gz"
    if !isfile(frequencyfile)
        println(frequencyfile, " does not exist.")
        continue
    end 
    freqs = parsefreqfile(frequencyfile)

    readsfile = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind$(truenumind)_basecount_IBS_ind$(truenumind).gz"
    if !isfile(readsfile)
        println(readsfile, " does not exist.")
        continue
    end 
    derivedreadsvec, coveragevec = parsereadfile(readsfile, ' ')
    #@show length(freqs) length(derivedreadsvec) length(coveragevec)
    @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
    
    t = tuple(derivedreadsvec, coveragevec, freqs)
    push!(results, t)
    println("Finished with anchor population ", anchorpop, " true number of individuals ", truenumind)
end # (anchorpop, ancientpop, truenumind)

for index in eachindex(allcombinations, results)
    println(allcombinations[index])
    reads, coverage, freq = results[index]
    println(rpad("Naive ", 25), length(reads))
    println(rpad("Without rounding ", 25), length(countmap(zip(reads, coverage, freq))))
    freq2 = round.(freq, digits=3)
    println("Rounding to three digits ", length(countmap(zip(reads, coverage, freq2))))
    freq2 = round.(freq, digits=2)
    println(rpad("Rounding to two digits ", 25), length(countmap(zip(reads, coverage, freq2))))
    freq2 = round.(freq, digits=1)
    println(rpad("Rounding to one digit ", 25), length(countmap(zip(reads, coverage, freq2))))
end 