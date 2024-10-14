# MCMCknowntau.jl
# Author: Jan van Waaij
# Created: 8 May 2024
# Last update: 8 May 2024
#
# Purpose:
# Run MCMC with the files data/frequencies/freqanchor/two_demes_demes_freqanchor.txt.gz 
# and 
# data/reads/two_demes_demes_1_2_reads_with_errors.txt.gz (1 individual, coverage 2), 
# number of individuals 1..8 and coverage 2..20. 
#
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia --optimize=3 path/to/MCMCknowntau.jl.
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

# nodes 01 en 10 doen het niet. 
nodes = [
    ("compute01", 24), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32
    ("compute02", 47), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32
    ("compute03", 29), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32 
]
sum(node[2] for node in nodes)

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/"
project = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/scripts/subsamplingMCMC/"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
numthreadsperprocess = 4
exeflags = ["--project=$project", "--threads=$numthreadsperprocess", "--optimize=3", "--banner=no", "--warn-overwrite=no"] 

addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

@everywhere using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed



nsteps = 100_000 # number of steps in the mcmc chain 
prioronn = DiscreteUniform(1,10)
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 

trueϵ = 0.003 # true error rate used in the simulations. 

for numind in 3:8
    coverages = 100:-1:1
    allcombinations = [(numind, cov) for cov in coverages]

    println("starting with ", numind, " individuals.")
    @sync @distributed for (numberofindividualsinthesediment, coverage) in allcombinations   
        local starttime = now()
        frequencyfile = "data/frequencies/freqanchor/two_demes_demes_freqanchor.txt.gz"
        if !isfile(frequencyfile)
            println(frequencyfile, " does not exist.")
            continue
        end 
        freqs = parsefreqfile(frequencyfile)
        readsfile = string("data/reads/subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", coverage, "_reads_with_errors.txt.gz")
        if !isfile(readsfile)
            println(readsfile, " does not exist.")
            continue
        end 
        derivedreadsvec, coveragevec = parsereadfile(readsfile)
        @assert length(unique(coveragevec)) == 1
        @assert coveragevec[1] == coverage
        @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

        csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", coverage, "_MCMC.csv"))
        zippedcsvfile = csvfile * ".gz"
        df = DataFrame(
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
        CSV.write(csvfile, df, names = [
            "trueϵ", 
            "truenumind", 
            "coverage", 
            "nsample", 
            "τCsample", 
            "τAsample", 
            "ϵsample", 
            "accepted", 
            "logjointprob", 
            "chainid"
        ])
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
                everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) 
                acceptedvecunpacked[index] = !everythingthesame
            end 
            local df = DataFrame(
                    trueϵ = fill(trueϵ, nsteps),
                    truenumind = fill(numberofindividualsinthesediment, nsteps),
                    coverage = fill(coverage, nsteps),
                    nsample = nsample,
                    τCsample = τCsampleunpacked,
                    τAsample = τAsampleunpacked,
                    ϵsample = ϵsampleunpacked,
                    accepted = acceptedvecunpacked,
                    logjointprob = logjointprob[:,1],
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
    println("finished with ", numind, " individuals.")
end 

finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("It took ", numhours, " and ", remainingtimeminutes, " to finish all files!")
