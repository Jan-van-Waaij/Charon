# runmcmconallfiles.jl
# Author: Jan van Waaij
# Created: 7 September 2023
# Last update: 11 September 2023
#
# Purpose:
# Run MCMC with errors.
#
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia -O3 --threads=auto path/to/runmcmconallfiles.jl
#
# Requirement:
# Frequency files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqgenotypes
# Read files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads
# 
# Output:
# ????
# The CSV file has 10 columns: NeAnchorPop, NeAncientPop, jointtimeyears, trueτC, trueτA, truenumind, τC, τA, numind, loglik.


using Data, # For parsing frequency and read files. 
      #GZip, # Read and write zip files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads # for threads

#=
Wat gaan we doen? 

* voor elke anchor gaan we tauA en tauC berekenen tussen de anchor en de ancient populatie. 
* We doen dat dan ook nog voor de verschillende coverages en number of individuals
=#

infoblock = 9
ancientpop = infoblock .+ (1:5000)
anchorpops = [infoblock .+ (1+5000:2*5000), infoblock .+ (1+2*5000:3*5000), infoblock .+ (1+3*5000:4*5000), infoblock .+ (1+4*5000:5*5000)]
trueτAs = (1000:-200:400) ./ (2*2000) 
trueτCs = trueτAs
jointtimegenerations = 1000:-200:400

nsteps = 100_000
prioronτCτA = Product([Uniform(), Uniform()])
prioronϵ = Uniform(0, 0.5)
nchains = 4

trueϵ = 0.003


@threads for maxcov in 2:20 # 7:4:20
    @threads for maxind in 2:8 
        starttime = now()
        numberofindividualsinthesedimentvec = Base.OneTo(maxind) #1:5 #  [1, 2, 4, 8]  #  [1, 2, 4, 8, 16, 32] #, 16, 32, 64] 
        coverage= maxcov#10    
        csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("multipledemestraditionalmcmccov_", maxcov, "_numind_", maxind, ".csv"))
        zippedcsvfile = csvfile * ".gz"
        if isfile(zippedcsvfile)
            println(zippedcsvfile, " does already exist. Do not overwrite existing files!")
            continue 
        end

        prioronn = DiscreteUniform(1, maximum(numberofindividualsinthesedimentvec))

        df = DataFrame(anchorpopulation = Int[],
                                NeAnchorPop = Int[],
                                NeAncientPop = Int[],
                                jointtimegenerations = Int[],
                                trueτC = Float64[],
                                trueτA = Float64[],
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
        for anchorpopulation in 1:4
            frequencyfile = "data/frequencies/freqanchor/multiple_demes_freq_anchor$(anchorpopulation).txt.gz"
            freqs = parsefreqfile(frequencyfile)
            
            trueτA = trueτAs[anchorpopulation]
            trueτC = trueτCs[anchorpopulation]
            for truenumind in numberofindividualsinthesedimentvec
                firstpart = "data/reads/multiple_demes_"
                readsfile = string(firstpart, 
                truenumind, '_', coverage, "_reads_with_errors.txt.gz")
                if !isfile(readsfile)
                    println(readsfile, " does not exist.")
                    continue
                end 
                derivedreadsvec, coveragevec = parsereadfile(readsfile)
                @assert length(unique(coveragevec)) == 1
                @assert coveragevec[1] == coverage
                @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
        # MCMCSimProp(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
                # NEW STUFF
                results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 100_000)
                for i in eachindex(results)
                    nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[i]
                    local df = DataFrame(anchorpopulation = repeat([anchorpopulation], nsteps),
                                NeAnchorPop = repeat([2000], nsteps),
                                NeAncientPop = repeat([2000], nsteps),
                                jointtimegenerations = repeat([jointtimegenerations[anchorpopulation]], nsteps),
                                trueτC = repeat([trueτC], nsteps),
                                trueτA = repeat([trueτA], nsteps),
                                trueϵ = repeat([trueϵ], nsteps),
                                truenumind = repeat([truenumind], nsteps),
                                coverage = repeat([coverage], nsteps),
                                nsample = nsample,
                                τCsample = τCsample,
                                τAsample = τAsample,
                                ϵsample = ϵsample,
                                accepted = acceptedvec,
                                logjointprob = logjointprob,
                                chainid = repeat([Int(i)], nsteps))
                    CSV.write(csvfile, df, append=true)
                end
                println("finished true num ind ", truenumind)
            end # truenumind 
        end # anchorpopulation
        command = `gzip $csvfile`
        run(command)
        finishtime = now()
        println(finishtime - starttime)
        println("Finished with ", csvfile)
    end # maxind  
end # maxcov