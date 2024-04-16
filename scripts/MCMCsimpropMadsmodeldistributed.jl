# MCMCsimpropMadsmodeldistributed.jl
# Author: Jan van Waaij
# Created: 12 April 2024
# Last update: 12 April 2024
#
# Purpose:
# Run MCMC with the files simulated by Mads.
#
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia --optimize=3 path/to/MCMCsimpropMadsmodeldistributed.jl
#

using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed

global starttime = now()

ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec] # 64 elementen. 
allcombinations = allcombinations[1:1]

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--threads=4", "--optimize=3"]

nodes = [
    ("compute01", 32), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32 
    #("node05", 6), #6), # 2 cores in use. MAX CAP FOR ME:5 ASO MODE 6
    ("node03", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    #("node04", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node12", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4 
    ("node01", 6), # not used. MAX CAP FOR ME: 5 ASO MODE 6
#    ("node02", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("node06", 4), # 2 core in use. MAX CAP FOR ME: 5 ASO MODE 6
    ("node07", 12), #14), # 16), # fully used MAX CAP FOR ME: 13 ASO MODE 16
    #("node13", 4), # 2 used. MAX CAP FOR ME: 3 ASO MODE 4
    #("node11", 10), #10), # fully used MAX CAP FOR ME: 8 ASO MODE 10
    #("node14", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    #("node08", 5), # not used MAX CAP FOR ME: 4 ASO MODE 5
    
]
sum(node[2] for node in nodes)

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
numthreadsperprocess = 4
exeflags = ["--project=$dir", "--threads=$numthreadsperprocess", "--optimize=3", "--banner=no"] 

#=
addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

@everywhere using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed
=#

nsteps = 100_000 # number of steps in mcmc chain 
prioronn = DiscreteUniform(1,10)
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 

trueϵ = 0.001 # true error rate used in the simulations. 

coverage = 64 # the coverage is 64 everywhere. 

#@sync @distributed 
for (anchorpop, ancientpop, truenumind) in allcombinations   
    local starttime = now()
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

    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("simprop_Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
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
    results = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 10_000)
    for chain in eachindex(results)
        nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[chain]
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
            everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) 
            acceptedvecunpacked[index] = !everythingthesame
        end 
        local df = DataFrame(
                            anchorpop = repeat([anchorpop], nsteps),
                            ancientpop = repeat([ancientpop], nsteps),
                            trueϵ = repeat([trueϵ], nsteps),
                            truenumind = repeat([truenumind], nsteps),
                            coverage = repeat([coverage], nsteps),
                            nsample = nsample,
                            τCsample = τCsampleunpacked,
                            τAsample = τAsampleunpacked,
                            ϵsample = ϵsampleunpacked,
                            accepted = acceptedvecunpacked,
                            logjointprob = logjointprob[:, 1],
                            chainid = fill(Int(chain), nsteps)
                            )
        CSV.write(csvfile, df, append=true)
    end # chain 
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
