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
      Base.Threads, # for threads
      Distributed

global starttime = now()

# /home/projects/gabriel/Jan_allelefreq/hg38/bams/IBS_ind16_basecount_IBS_ind16.gz
# let op de coverage is niet constant. 
anchorvec = 1:4
numindvec = 1:8
coveragevec = [10, 15, 20] #2:20
allcombinations = [(anchor, numind, coverage) for anchor in anchorvec for numind in numindvec for coverage in coveragevec] # 608 num ind 

nsteps = 100_000 # number of steps in mcmc chain 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 

trueϵ = 0.001 # true error rate used in the simulations. 

nodes = [
    ("compute01", 20), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32 
    #("node05", 6), #6), # 2 cores in use. MAX CAP FOR ME:5 ASO MODE 6
    #("node03", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    #("node04", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node12", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4 
    ("node01", 6), # not used. MAX CAP FOR ME: 5 ASO MODE 6
    ("node02", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("node06", 6), # 2 core in use. MAX CAP FOR ME: 5 ASO MODE 6
    #("node07", 16), #14), # 16), # fully used MAX CAP FOR ME: 13 ASO MODE 16
    ("node13", 4), # 2 used. MAX CAP FOR ME: 3 ASO MODE 4
    #("node11", 10), #10), # fully used MAX CAP FOR ME: 8 ASO MODE 10
    ("node14", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("node08", 5), # not used MAX CAP FOR ME: 4 ASO MODE 5
]
#nodes = [("compute01", 1)]
sum(node[2] for node in nodes)

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
numthreadsperprocess = 4
exeflags = ["--project=$dir", "--threads=$numthreadsperprocess", "--optimize=3", "--banner=no"] 

addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

@everywhere using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Distributed # for distributed computing.

@sync @distributed for (anchor, numind, coverage) in allcombinations  
    local starttime = now()
    frequencyfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqanchor/multiple_demes_freq_anchor$anchor.txt.gz"
    if !isfile(frequencyfile)
        println(frequencyfile, " does not exist.")
        continue
    end 
    freqs = parsefreqfile(frequencyfile)
    readsfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads/multiple_demes_$(numind)_$(coverage)_reads_with_errors.txt.gz"
    if !isfile(readsfile)
        println(readsfile, " does not exist.")
        continue
    end 
    derivedreadsvec, coveragevec = parsereadfile(readsfile)
    @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
    @assert length(unique(coveragevec)) == 1 
    @assert coveragevec[1] == coverage
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "mcmc_multiple_demes_anchor_$(anchor)_num_ind_$(numind)_cov_$(coverage).csv")
    zippedcsvfile = csvfile * ".gz"
    #=if isfile(zippedcsvfile)
        println(zippedcsvfile, " does already exist. Do not overwrite existing files!")
        continue 
    end=#
    df = DataFrame(anchorpop = Int[],
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
    results = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 0, false)
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
            everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) && (nsample[index] == nsample[index-1])
            acceptedvecunpacked[index] = !everythingthesame
        end 
        local df = DataFrame(
                            anchorpop = fill(anchor, nsteps),
                            trueϵ = fill(trueϵ, nsteps),
                            truenumind = fill(numind, nsteps),
                            coverage = fill(coverage, nsteps),
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
    try 
        command = `gzip -f $csvfile`
        run(command)
    catch
        println("I was not able to zip $csvfile.")
    end 
    println("I am finished with ", csvfile)
    finishtime = now()
    timedifference = finishtime - starttime
    numhours = round(timedifference, Hour, RoundDown)
    remainingtime = timedifference - numhours
    remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
    remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
    println("It took me ", numhours, " and ", remainingtimeminutes, '.')
end # (anchorpop, ancientpop, truenumind)

finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("It took me ", numhours, " and ", remainingtimeminutes, " to finish all files!")

rmprocs(workers())

