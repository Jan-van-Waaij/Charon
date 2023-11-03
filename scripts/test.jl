# runmcmconallfiles.jl
# Author: Jan van Waaij
# Created: 4 May 2023
# Last update: 4 May 2023
#
# Purpose:
# Run MCMC.
#
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia path/to/runmcmconallfiles.jl
#
# Requirement:
# Frequency files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqgenotypes
# Read files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads
# 
# Output:
# ????
# The CSV file has 10 columns: NeAnchorPop, NeAncientPop, jointtimeyears, trueτC, trueτA, truenumind, τC, τA, numind, loglik.

#using Pkg
#Pkg.update()

using Data, # For parsing frequency and read files. 
      GZip, # Read and write zip files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling # To sample from the posterior.

vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles")
allfiles = ["two_demes_3000_6000_30000_0.vcf.gz"]



fileswithoutextensionvec = [file[begin:end-7] for file in allfiles]

numberofindividualsinthesedimentvec =  [2] 

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "test.csv.gz")
# if isfile(csvfile)
#     println(csvfile, " does already exist. Do not overwrite existing files!")
#     exit(1)
# end
openedcsvfile = GZip.open(csvfile, "a") # the a of append. 

df = DataFrame(NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    trueτC = Float64[],
                    trueτA = Float64[],
                    truenumind = Int[],
                    τCsample = Float64[],
                    τAsample = Float64[],
                    accepted = Bool[]
                    )

#CSV.write(openedcsvfile, df)

nsteps = 100_000 
prior = Product([Uniform(), Uniform()])

#vcffilewithoutextension = "two_demes_11000_6000_70000_0"

# sp = split(vcffilewithoutextension, '_')
# NeAncientPop = parse(Int, sp[3])
# NeAnchorPop = parse(Int, sp[4])
# jointtimeyears = parse(Int, sp[5])
# generationtime = 25
# trueτA = jointtimeyears / (2*NeAncientPop*25)
# trueτC = jointtimeyears / (2*NeAnchorPop*25)
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
                string(filewithoutextension, '_', truenumind, "_100_reads.txt.gz") )
            if !isfile(readsfile)
                println(readsfile, " does not exist.")
                continue
            end 
            gzopen(readsfile) do f 
                global derivedreadsvec, coveragevec = parsereadfile(f)
            end
            @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
            
            init = (τC = rand(), τA = rand())

            @code_warntype MHsampler(nsteps, init, prior, coveragevec, derivedreadsvec, freqs, truenumind, 25_000) 

            #=global df = DataFrame(NeAnchorPop = repeat([NeAnchorPop], nsteps),
                        NeAncientPop = repeat([NeAncientPop], nsteps),
                        jointtimeyears = repeat([jointtimeyears], nsteps),
                        trueτC = repeat([trueτC], nsteps),
                        trueτA = repeat([trueτA], nsteps),
                        truenumind = repeat([truenumind], nsteps),
                        τCsample = τCsample,
                        τAsample = τAsample,
                        accepted = acceptedvec)
            CSV.write(openedcsvfile, df, append=true)
            =#
            println("finished ", filewithoutextension, " true num ind ", truenumind)
        end 
    end
finally
    close(openedcsvfile)
end 
