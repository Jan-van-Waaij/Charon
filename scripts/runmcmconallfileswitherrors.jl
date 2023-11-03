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

numberofindividualsinthesedimentvec = 1:3 #  [1, 2, 4, 8]  #  [1, 2, 4, 8, 16, 32] #, 16, 32, 64] 
maxreadsvec = [9]

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "allsampleswitherrorsthreepromileandcoverage9max3individuals.csv.gz")

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
prioronn = DiscreteUniform(1, 3) # Discrete uniform prior on 1,...,10.
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
            for maxreads in maxreadsvec
                readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", 
                    string(filewithoutextension, '_', truenumind, '_', maxreads, "_reads_with_errors.txt.gz"))
                if !isfile(readsfile)
                    println(readsfile, " does not exist.")
                    continue
                end 
                gzopen(readsfile) do f 
                    global derivedreadsvec, coveragevec = parsereadfile(f)
                end
                @assert length(unique(coveragevec)) == 1
                coverage = coveragevec[1]
                @assert coverage == maxreads
                @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
                results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs)
                for i in eachindex(results)
                    nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[i]
                    local df = DataFrame(NeAnchorPop = repeat([NeAnchorPop], nsteps),
                                NeAncientPop = repeat([NeAncientPop], nsteps),
                                jointtimeyears = repeat([jointtimeyears], nsteps),
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
                println("finished ", filewithoutextension, " true num ind ", truenumind)
            end # maxreads
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