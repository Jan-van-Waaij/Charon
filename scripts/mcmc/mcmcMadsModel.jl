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

using Distributed 

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
project = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/scripts/mcmc/"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$project", "--threads=4", "--optimize=3"]

nodes = [
	("compute01", 32),
	("compute02", 48),
	#("compute03", 10),
	("node02", 4),
	("node03", 6),
	("node04", 6),
	("node05", 6),
	#("node06", 6),
	#("node07", 4),
	("node08", 5),
	("node09", 6),
	("node10", 4),
	("node11", 10),
	("node12", 4),
	("node13", 4),
	("node14", 4)
]

addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

# trigger precompilation.
@everywhere begin 
    using Data # For parsing frequency and read files. 
    using DataFrames # work with data frames.
    using CSV # convert dataframes to a CSV file. 
    using Distributions # To specify the prior.
    using PosteriorSampling # To sample from the posterior.
    using Dates # for the now function 
    using Base.Threads # for threads
    using Distributed 
    using SharedArrays # for sharing an array among different processes. 
end 

starttime = now()

ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops =  ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = 1:8
coveragevec = setdiff(5:5:100, [100, 50, 20, 10])
allcombinations = [(anchorpop, ancientpop, truenumind, coverage) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec for coverage in coveragevec ] # 64 elementen. 
filter!(tuple -> !(tuple[1] == "CEU" && tuple[2] == "CEUS"), allcombinations)

nsteps = 100_000 # number of steps in mcmc chain 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 


trueϵ = 0.003 # true error rate used in the simulations. 

@sync @distributed for index in eachindex(allcombinations)
    anchorpop, ancientpop, truenumind, coverage = allcombinations[index]
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "MadsModel", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind,  "_cov_", coverage, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if isfile(zippedcsvfile)
        println("skipping $zippedcsvfile")
        continue
    end 

    starttime = now()
    frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", "freqanchor", "Madsmodel_freq_anchor_$(anchorpop).txt.gz")
    if !isfile(frequencyfile)
        println(frequencyfile, " does not exist.")
        continue
    end 
    freqs = parsefreqfile(frequencyfile)
    readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", string("Mads_model_pop_", ancientpop, "_num_ind_",
    truenumind, "_cov_", coverage, "_reads_with_errors.txt.gz"))
    if !isfile(readsfile)
        println(readsfile, " does not exist.")
        continue
    end 
    derivedreadsvec, coveragevec = parsereadfile(readsfile)
    @assert length(unique(coveragevec)) == 1
    @assert coveragevec[1] == coverage
    @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

    #=
    if isfile(zippedcsvfile)
        println(zippedcsvfile, " does already exist. Do not overwrite existing files!")
        continue 
    end
    =#
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
    #results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
    chains = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 0, false)
    for index in eachindex(chains)
        nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = chains[index]
        τCsampleunpacked = Vector{Float64}(undef, nsteps)
        τAsampleunpacked = Vector{Float64}(undef, nsteps)
        ϵsampleunpacked = Vector{Float64}(undef, nsteps)
        for index in 1:nsteps
            τCsampleunpacked[index] = τCsample[index, nsample[index]]
            τAsampleunpacked[index] = τAsample[index, nsample[index]]
            ϵsampleunpacked[index] = ϵsample[index, nsample[index]]
        end 
        acceptedvecunpacked = BitVector(undef, nsteps)
        acceptedvecunpacked[1] = true 
        for index in 2:nsteps
            everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) && (nsample[index] == nsample[index-1])
            acceptedvecunpacked[index] = !everythingthesame
        end 
        local df = DataFrame(
                            anchorpop = fill(anchorpop, nsteps),
                            ancientpop = fill(ancientpop, nsteps),
                            trueϵ = fill(trueϵ, nsteps),
                            truenumind = fill(truenumind, nsteps),
                            coverage = fill(coverage, nsteps),
                            nsample = nsample,
                            τCsample = τCsampleunpacked,
                            τAsample = τAsampleunpacked,
                            ϵsample = ϵsampleunpacked,
                            accepted = acceptedvecunpacked,
                            logjointprob = logjointprob[:,1],
                            chainid = fill(Int(index), nsteps)
                            )
        CSV.write(csvfile, df, append=true)
    end  
    command = `gzip -f $csvfile`
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

finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("It took ", numhours, " and ", remainingtimeminutes, " to finish all files!")

#=
Finished duurde ongeveer 30-38 uur. 
=#