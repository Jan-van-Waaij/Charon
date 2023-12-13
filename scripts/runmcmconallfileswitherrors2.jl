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
      GZip, # Read and write zip files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads # for threads


@threads for maxcov in 7:4:20 # 7:4:20
    @threads for maxind in 3:maxcov 
        starttime = now()
        numberofindividualsinthesedimentvec = Base.OneTo(maxind) #1:5 #  [1, 2, 4, 8]  #  [1, 2, 4, 8, 16, 32] #, 16, 32, 64] 
        coverage= maxcov#10
        frequencyfile = "data/frequencies/freqanchor/two_demes_demes_freqanchor.txt.gz"
        freqs = parsefreqfile(frequencyfile)

        csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("demesmcmccov", maxcov, "numind", maxind, ".csv.gz"))

        if isfile(csvfile)
            println(csvfile, " does already exist. Do not overwrite existing files!")
            continue 
        end

        df = DataFrame(NeAnchorPop = Int[],
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

        nsteps = 400_000 
        prioronτCτA = Product([Uniform(), Uniform()])
        prioronn = DiscreteUniform(1, maximum(numberofindividualsinthesedimentvec)) # Discrete uniform prior on 1,...,10.
        prioronϵ = Uniform(0, 0.5)
        nchains = 4

        #vcffilewithoutextension = "two_demes_11000_6000_70000_0"

        # sp = split(vcffilewithoutextension, '_')
        # NeAncientPop = parse(Int, sp[3])
        # NeAnchorPop = parse(Int, sp[4])
        # jointtimeyears = parse(Int, sp[5])
        # generationtime = 25
        # trueτA = jointtimeyears / (2*NeAncientPop*25)
        # trueτC = jointtimeyears / (2*NeAnchorPop*25)
        trueϵ = 0.003
        trueτA = 1000 / (2*2000)
        trueτC = 1000 / (2*2000)
        try
            for truenumind in numberofindividualsinthesedimentvec
                firstpart = "data/reads/two_demes_demes_"
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
                results = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs)
        #       results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs)
                for i in eachindex(results)
                    nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[i]
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
                    local df = DataFrame(NeAnchorPop = repeat([2000], nsteps),
                                NeAncientPop = repeat([1000], nsteps),
                                jointtimegenerations = repeat([1000], nsteps),
                                trueτC = repeat([trueτC], nsteps),
                                trueτA = repeat([trueτA], nsteps),
                                trueϵ = repeat([trueϵ], nsteps),
                                truenumind = repeat([truenumind], nsteps),
                                coverage = repeat([coverage], nsteps),
                                nsample = nsample,
                                τCsample = τCsampleunpacked,
                                τAsample = τAsampleunpacked,
                                ϵsample = ϵsampleunpacked,
                                accepted = acceptedvecunpacked,
                                logjointprob = logjointprob[:,1],
                                chainid = repeat([Int(i)], nsteps))
                    CSV.write(csvfile, df, append=true)
                end 
                println("finished true num ind ", truenumind)
            end # truenumind
        finally
            #close(openedcsvfile)
            finishtime = now()
            println(finishtime - starttime)
            println("Finished with ", csvfile)
        end # try 
    end 
end 
# 