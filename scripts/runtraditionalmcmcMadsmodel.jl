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
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads # for threads

ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec]

nsteps = 200_000 # number of steps in mcmc chain 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 

trueϵ = 0.001 # true error rate used in the simulations. 

coverage = 64 # the coverage is 64 everywhere. 

@threads for (anchorpop, ancientpop, truenumind) in allcombinations    
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

    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if isfile(zippedcsvfile)
        println(zippedcsvfile, " does already exist. Do not overwrite existing files!")
        continue 
    end
    df = DataFrame(anchorpop = String[],
                   ancientpop = String[],
                   trueϵ = Float64[],
                   truenumind = Int[],
                   coverage = Int[],
                   nsample = Int[],
                   τCsample = Float64[],
                   τAsample = Float64[],
                   ϵsample = Float64[],
                   accepted = Bool[],
                   logjointprob = Float64[],
                   chainid = Int[]
                  )
    CSV.write(csvfile, df)
    results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
    for chain in eachindex(results)
        nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[chain]
        #= 
        anchorpop = String[],
        ancientpop = String[],
        trueϵ = Float64[],
        truenumind = Int[],
        coverage = Int[],
        nsample = Int[],
        τCsample = Float64[],
        τAsample = Float64[],
        ϵsample = Float64[],
        accepted = Bool[],
        logjointprob = Float64[],
        chainid = Int[]
        =#
        local df = DataFrame(
                             anchorpop = repeat([anchorpop], nsteps),
                             ancientpop = repeat([ancientpop], nsteps),
                             trueϵ = repeat([trueϵ], nsteps),
                             truenumind = repeat([truenumind], nsteps),
                             coverage = repeat([coverage], nsteps),
                             nsample = nsample,
                             τCsample = τCsample,
                             τAsample = τAsample,
                             ϵsample = ϵsample,
                             accepted = acceptedvec,
                             logjointprob = logjointprob,
                             chainid = repeat([Int(chain)], nsteps)
                            )
        CSV.write(csvfile, df, append=true)
    end # chain 
    command = `gzip $csvfile`
    run(command)
    finishtime = now()
    timedifference = finishtime - starttime
    numhours = round(timedifference, Hour, RoundDown)
    remainingtime = timedifference - numhours
    remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
    remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
    println("It took ", numhours, " and ", remainingtimeminutes, '.')
    println("Finished with ", csvfile)
end # (anchorpop, ancientpop, truenumind)

#=
Finished duurde ongeveer 30-38 uur. 
=#

#=
@threads for anchorpop in anchorpops    
    frequencyfile = joinpath("/", "home", "projects", "gabriel", "Jan_allelefreq", "hg38", "bams", string("IBS_ind1_freq_", anchorpop, ".gz"))
    if !isfile(frequencyfile)
        println(frequencyfile, " does not exist.")
        continue
    end 
    freqs = parsefreqfile(frequencyfile)
    @threads for ancientpop in ancientpops    
        @threads for truenumind in numberofindividualsinthesedimentvec
            starttime = now()
            readsfile = joinpath("/", "home", "projects", "gabriel", "Jan_allelefreq", "simulations", "split", string("all_basecount_", ancientpop, '_', truenumind, ".gz"))
            if !isfile(readsfile)
                println(readsfile, " does not exist.")
                continue
            end 
            derivedreadsvec, coveragevec = parsereadfile(readsfile, ' ')
            @assert length(unique(coveragevec)) == 1
            @assert coveragevec[1] == coverage
            @show length(freqs) length(derivedreadsvec) length(coveragevec)
            @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

            csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
            zippedcsvfile = csvfile * ".gz"
            if isfile(zippedcsvfile)
                println(zippedcsvfile, " does already exist. Do not overwrite existing files!")
                continue 
            end
            df = DataFrame(anchorpop = String[],
                           ancientpop = String[],
                           trueϵ = Float64[],
                           truenumind = Int[],
                           coverage = Int[],
                           nsample = Int[],
                           τCsample = Float64[],
                           τAsample = Float64[],
                           ϵsample = Float64[],
                           accepted = Bool[],
                           logjointprob = Float64[],
                           chainid = Int[]
                          )
            CSV.write(csvfile, df)
            results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 0)
            for chain in eachindex(results)
                nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[chain]
                #= 
                anchorpop = String[],
                ancientpop = String[],
                trueϵ = Float64[],
                truenumind = Int[],
                coverage = Int[],
                nsample = Int[],
                τCsample = Float64[],
                τAsample = Float64[],
                ϵsample = Float64[],
                accepted = Bool[],
                logjointprob = Float64[],
                chainid = Int[]
                =#
                local df = DataFrame(anchorpop = repeat([anchorpop], nsteps),
                                     ancientpop = repeat([ancientpop], nsteps),
                                     trueϵ = repeat([trueϵ], nsteps),
                                     truenumind = repeat([truenumind], nsteps),
                                     coverage = repeat([coverage], nsteps),
                                     nsample = nsample,
                                     τCsample = τCsample,
                                     τAsample = τAsample,
                                     ϵsample = ϵsample,
                                     accepted = acceptedvec,
                                     logjointprob = logjointprob,
                                     chainid = repeat([Int(chain)], nsteps))
                CSV.write(csvfile, df, append=true)
            end # chain 
            command = `gzip $csvfile`
            run(command)
            finishtime = now()
            println(finishtime - starttime)
            println("Finished with ", csvfile)
        end # truenumind
    end # ancientpop
end # anchorpop 
=#