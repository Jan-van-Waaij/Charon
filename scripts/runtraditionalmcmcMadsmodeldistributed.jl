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
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed


ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec] # 64 elementen. 

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--threads=4", "--optimize=3"]

#= max nodes 
nodes = [
    ("node01", 6),
    ("node02", 4),
    ("node03", 6),
    ("node04", 6),
    ("node05", 6),
    ("node06", 6),
    ("node07", 16),
    ("node08", 5),
    ("node11", 10),
    ("node12", 4),
    ("node13", 4),
    ("node14", 4),
]
=#

nodes = [
    ("node01", 6), # not used
    ("node02", 4), # not used
    ("node03", 6), # not used 
    ("node04", 6), # not used 
    ("node05", 3), #6), # 3 cores in use. 
    ("node06", 6), # not used 
    ("node07", 13), #14), # 16), # 2 cores in use. 
    ("node08", 5), # not used 
    ("node11", 7), #10), # 3 cores in use 
    ("node12", 4), # not used 
    #("node13", 4), # 100 % in use 
    ("node14", 4), # not used 
]

numprocs = sum(item[2] for item in nodes)

@assert numprocs == 64

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
numthreadsperprocess = 4
exeflags = ["--project=$dir", "--threads=$numthreadsperprocess", "--optimize=3", "--banner=no"] 

processes = Int[]
for node in nodes 
    append!(processes, addprocs([node], dir=dir, exename=exename, exeflags=exeflags))
end 

@everywhere using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed

@assert nworkers() == 64
@assert nworkers() == length(processes)
@assert numprocs == 64
futures = [@spawnat process Threads.nthreads() for process in processes]
num_threads = fetch.(futures)
#@show num_threads # weghalen en nog een test toevoegen om te bepalen of de project goed is. 
@assert unique(num_threads) == [numthreadsperprocess]
futures2 = [@spawnat process Base.active_project() for process in processes]
projects = fetch.(futures2)
@assert unique(projects) == [dir*"/Project.toml"]


nsteps = 200_000 # number of steps in mcmc chain 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 
@assert nchains == numthreadsperprocess


trueϵ = 0.001 # true error rate used in the simulations. 

coverage = 64 # the coverage is 64 everywhere. 

@sync @distributed for (anchorpop, ancientpop, truenumind) in allcombinations    
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

#=
Finished duurde ongeveer 30-38 uur. 
=#