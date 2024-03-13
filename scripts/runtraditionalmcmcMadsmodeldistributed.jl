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

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--threads=4", "--optimize=3"]

nodes = [
    ("node01", 6), # not used. MAX CAP FOR ME: 5 ASO MODE 6
    ("node02", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("node03", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node04", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node05", 5), #6), # 2 cores in use. MAX CAP FOR ME:5 ASO MODE 6
    ("node06", 5), # 1 core in use. MAX CAP FOR ME: 5 ASO MODE 6
    #("node07", 0), #14), # 16), # fully used MAX CAP FOR ME: 13 ASO MODE 16
    ("node08", 5), # not used MAX CAP FOR ME: 4 ASO MODE 5
    #("node11", 0), #10), # fully used MAX CAP FOR ME: 8 ASO MODE 10
    ("node12", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4 
    ("node13", 3), # 2 used. MAX CAP FOR ME: 3 ASO MODE 4
    ("node14", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("compute01", 16) # 50 percent used. MAX CAP FOR ME: 26 ASO MODE 32 
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

@everywhere function forceprecompile()
    nchains = 2
    nsteps = 4 
    prioronn = DiscreteUniform(1, 2)
    prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
    prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
    coverage = 4 
    derivedreadsvec = [1, 1, 3, 2]
    freqs = [0.1, 0.2, 0.1, 0.15]
    MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 1)
    return nothing 
end 

forceprecompile()
    

futures = [@spawnat process forceprecompile() for process in processes]
fetch.(futures)

@assert nworkers() == 64
@assert nworkers() == length(processes)
@assert numprocs == 64
futures2 = [@spawnat process Threads.nthreads() for process in processes]
num_threads = fetch.(futures2)
#@show num_threads # weghalen en nog een test toevoegen om te bepalen of de project goed is. 
@assert unique(num_threads) == [numthreadsperprocess]
futures3 = [@spawnat process Base.active_project() for process in processes]
projects = fetch.(futures3)
@assert unique(projects) == [dir*"/Project.toml"]


nsteps = 400_000 # number of steps in mcmc chain 
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
    CSV.write(csvfile, df, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])
    #results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
    results = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
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
            if (nsample[index] == nsample[index-1]) && (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1])  
                acceptedvecunpacked[index] = false
            else
                acceptedvecunpacked[index] = true
            end 
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
                             logjointprob = logjointprob[:,1],
                             chainid = repeat([Int(chain)], nsteps)
                            )
        CSV.write(csvfile, df, append=true, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])
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