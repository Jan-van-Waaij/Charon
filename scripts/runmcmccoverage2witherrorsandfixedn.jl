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
      Dates # for the now function 

starttime = now()
vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles")
allfilesunfiltered = readdir(vcffilespath)
allfiles = String[]
for file in allfilesunfiltered 
    if file[end-6:end] != ".vcf.gz" # remove non VCF files.
        continue
    end 
    sp = split(file, '_')
    if sp[3] == "1000" # remove two_demes_1000_10_0.vcf.gz
        continue 
    end 
    push!(allfiles, file)
end 

fileswithoutextensionvec = [file[begin:end-7] for file in allfiles]

numberofindividualsinthesedimentvec = 1:8 #  [1, 2, 4, 8]  #  [1, 2, 4, 8, 16, 32] #, 16, 32, 64] 

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "allsamplescoverage2fixednwitherrors.csv.gz")

if isfile(csvfile)
    println(csvfile, " does already exist. Do not overwrite existing files!")
    exit(1)
end

#openedcsvfile = GZip.open(csvfile, "a") # the a of append. 

df = DataFrame(NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    trueτC = Float64[],
                    trueτA = Float64[],
                    trueϵ = Float64[],
                    truenumind = Int[],
                    τCsample = Float64[],
                    τAsample = Float64[],
                    ϵsample = Float64[],
                    accepted = Bool[],
                    logjointprob = Float64[],
                    chainid = Int[]
                    )

#CSV.write(openedcsvfile, df)
CSV.write(csvfile, df)

nsteps = 100_000 
prioronτCτA = Product([Uniform(), Uniform()])
#prioronn = DiscreteUniform(1, 2) # Discrete uniform prior on 1, 2, 3.
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
try
    for filewithoutextension in fileswithoutextensionvec
        # filewithoutextension is of the form two_demes_Ne1_Ne2_jointimeyears_0
        sp = split(filewithoutextension, '_')
        NeAncientPop = parse(Int, sp[3])
        NeAnchorPop = parse(Int, sp[4])
        jointtimeyears = parse(Int, sp[5])
        generationtime = 25
        trueτA = jointtimeyears / (2*NeAncientPop*generationtime)
        trueτC = jointtimeyears / (2*NeAnchorPop*generationtime)
        frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", 
        "freqanchor", string(filewithoutextension, "_freqanchor.txt.gz"))
        gzopen(frequencyfile) do f 
            global freqs = parsefreqfile(f)
        end
        for truenumind in numberofindividualsinthesedimentvec
            readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", 
                string(filewithoutextension, '_', truenumind, '_', "2_reads_with_errors.txt.gz"))
            if !isfile(readsfile)
                println(readsfile, " does not exist.")
                continue
            end 
            gzopen(readsfile) do f 
                global derivedreadsvec, coveragevec = parsereadfile(f)
            end
            @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
            @assert unique(coveragevec) == [2]
            coverage = 2
# MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
            results = MHsampler(nchains, nsteps, truenumind, prioronτCτA, prioronϵ, coveragevec , derivedreadsvec, freqs)
            #results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, coverage, derivedreadsvec, freqs) 
            for i in eachindex(results)
# (nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob)
                τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[i]
# (NeAnchorPop = Int[],
# NeAncientPop = Int[],
# jointtimeyears = Int[],
# trueτC = Float64[],
# trueτA = Float64[],
# trueϵ = Float64[],
# truenumind = Int[],
# coverage = Int[],
# nsample = Int[],
# τCsample = Float64[],
# τAsample = Float64[],
# ϵsample = Float64[],
# accepted = Bool[],
# logjointprob = Float64[],
# chainid = Int[]
# )
                local df = DataFrame(NeAnchorPop = repeat([NeAnchorPop], nsteps),
                            NeAncientPop = repeat([NeAncientPop], nsteps),
                            jointtimeyears = repeat([jointtimeyears], nsteps),
                            trueτC = repeat([trueτC], nsteps),
                            trueτA = repeat([trueτA], nsteps),
                            trueϵ = repeat([trueϵ], nsteps),
                            truenumind = repeat([truenumind], nsteps),
                            τCsample = τCsample,
                            τAsample = τAsample,
                            ϵsample = ϵsample,
                            accepted = acceptedvec,
                            logjointprob = logjointprob,
                            chainid = repeat([Int(i)], nsteps))
                #CSV.write(openedcsvfile, df, append=true)
                CSV.write(csvfile, df, append=true)
            end # i in results.
            println("finished ", filewithoutextension, " true num ind ", truenumind)
        end # truenumind
    end # filewithoutextension
finally
    #close(openedcsvfile)
    finishtime = now()
    println(finishtime - starttime)
    println("Finished with ", csvfile)
end # try 

# First run. 22.4 uur. ERGUH LANG.
# Initialising is not the problem. 
# Nieuwe versie met losse arrays. Nu 1.12 uur. 
# Multiple chains version: 1.43 uur.
# 6 uur en 24 minuten.