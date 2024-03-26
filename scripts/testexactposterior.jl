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

# trigger precompilation.
using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates # for the now function 
     
ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec] # 64 elementen. 

#randomselection = sample(1:64, 4, replace=false)

#= 
[38, 63, 36, 57]
("CHB", "YRIS", 2)
("KAR", "KARS", 4)
("CHB", "CHBS", 8)
("KAR", "CEUS", 1)

println(randomselection)
allcombinations = allcombinations[randomselection]
for combination in allcombinations
    println(combination)
end 
=# 


trueϵ = 0.001 # true error rate used in the simulations. 

coverage = 64 # the coverage is 64 everywhere. 

τrange = 0.1:0.1:1.0

allcombinations2 = allcombinations[[38]]
for (anchorpop, ancientpop, truenumind) in allcombinations2
    starttime = now()
    frequencyfile = "data/frequencies/freqanchor/Madsmodel_freq_anchor_$(anchorpop).txt.gz"
    if !isfile(frequencyfile)
        println(frequencyfile, " does not exist.")
        continue
    end 
    freqs = parsefreqfile(frequencyfile)
    readsfile = string("data/reads/Mads_model_pop_", ancientpop, "_num_ind_",
    truenumind, "_cov_", coverage, "_reads_with_errors.txt.gz")
    if !isfile(readsfile)
        println(readsfile, " does not exist.")
        continue
    end 
    derivedreadsvec, coveragevec = parsereadfile(readsfile)
    @assert length(unique(coveragevec)) == 1
    @assert coveragevec[1] == coverage
    @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

    results = ExactPosterior(truenumind, τrange, 0.001, coverage, derivedreadsvec, freqs, 1)
    df = DataFrame(τC = results[1,:], τA=results[2,:], logpost = results[3,:])
    println(df)
    csvfile = "test.csv"
    CSV.write(csvfile, df)
    command = `gzip -f $csvfile`
    run(command)
end # (anchorpop, ancientpop, truenumind)

