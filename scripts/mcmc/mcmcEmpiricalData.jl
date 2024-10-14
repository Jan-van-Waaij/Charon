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

# nice -19 /home/ctools/julia-1.6.7/bin/julia --project=scripts/mcmc --optimize=3 scripts/mcmc/mcmcEmpiricalData.jl

using Distributed 

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
project = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/scripts/mcmc/"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$project", "--threads=4", "--optimize=3"]

nodes = [
	("compute01", 32),
	("compute02", 36),
	("compute03", 16),
	("node02", 4),
	("node03", 6),
	#("node04", 6),
	("node05", 3),
	("node06", 6),
	("node07", 12),
	("node08", 5),
	("node09", 6),
	("node10", 4),
	("node11", 6),
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
end 

starttime = now()

folder = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/"
files = readdir(folder)
freqfiles = files[contains.(files, r"IBS_ind\d_freq_[A-Z]{3}.gz")]
populations = unique([freqfile[end-5:end-3] for freqfile in freqfiles])
allthecombinations = [(population, numind) for population in populations for numind in [1, 2, 4, 8]]
allthecombinations = allthecombinations[2:end]

nsteps = 100_000 # number of steps in mcmc chain 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 

csvfilefolder = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "EmpiricalData")
    

@sync @distributed for index in eachindex(allthecombinations)
    starttime = now()
    
    population, numind = allthecombinations[index]
    
    readfile = folder * "IBS_ind$(numind)_basecount_IBS_ind$(numind).gz"
    freqfile = folder * "IBS_ind$(numind)_freq_$(population).gz"
    csvfile = joinpath(csvfilefolder, string("EmpiricalData", "_anchor_pop_", population, "_true_num_ind_", numind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    
    if isfile(zippedcsvfile)
        println("skipping $zippedcsvfile")
        continue
    end 
    if !isfile(freqfile)
        println(freqfile, " does not exist.")
        continue
    end
    if !isfile(readfile)
        println(readfile, " does not exist.")
        continue
    end 

    freqs = parsefreqfile(freqfile)
    derivedreadsvec, coveragevec = parsereadfile(readfile, ' ')

    @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

    df = DataFrame(anchorpop = String[],
                truenumind = Int[],
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
    chains = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coveragevec, derivedreadsvec, freqs, 0, false)
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
                            anchorpop = fill(population, nsteps),
                            truenumind = fill(numind, nsteps),
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
    files = readdir(csvfilefolder)
    gzfiles = endswith.(files, ".csv.gz")
    println(round(Int, sum(gzfiles)/length(allthecombinations)*100), "% of files finished!")
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