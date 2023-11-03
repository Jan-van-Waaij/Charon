# runmcmconallfiles.jl
# Author: Jan van Waaij
# Created: 23 June 2023
# Last update: 23 June 2023
#
# Purpose:
# Run MCMC.
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

numberofindividualsinthesedimentvec = [1, 2, 4, 8]  #  [1, 2, 4, 8, 16, 32] #, 16, 32, 64] 
maxreadsvec = [5, 10, 15, 20, 100]

nrange = 1:10 
τCrange = 0.01:0.01:0.99
τArange = 0.01:0.01:0.99

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "maxlik", "maxlikoverallfiles.csv.gz")

if isfile(csvfile)
    println(csvfile, " does already exist. Do not overwrite existing files!")
    exit(1)
end

openedcsvfile = GZip.open(csvfile, "w") # the w write. 

df = DataFrame(NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    coverage = Int[],
                    trueτC = Float64[],
                    maxlikτC = Float64[],
                    trueτA = Float64[],
                    maxlikτA = Float64[],
                    truenumind = Int[],
                    maxliknumind = Int[]
                    )

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
                    string(filewithoutextension, '_', truenumind, '_', maxreads, "_reads.txt.gz"))
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

                # HIER KOMT DE MAXIMISER
                maxliknumind, maxlikτC, maxlikτA = maxlik(nrange, τCrange, τArange, coverage, derivedreadsvec, freqs)
                # push to the dataframe 
                push!(df, (NeAnchorPop,
                           NeAncientPop, 
                           jointtimeyears,
                           coverage,
                           trueτC,
                           maxlikτC,
                           trueτA,
                           maxlikτA,
                           truenumind,
                           maxliknumind)
                     )
#=
NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    coverage = Int[],
                    trueτC = Float64[],
                    maxlikτC = Float64[],
                    trueτA = Float64[],
                    maxlikτA = Float64[],
                    truenumind = Int[],
                    maxliknumind = Int[]
=#
                println("finished ", filewithoutextension, " true num ind ", truenumind)
            end # maxreads
        end # truenumind
    end # filewithoutextension
finally
    CSV.write(openedcsvfile, df)
    close(openedcsvfile)
    finishtime = now()
    println(finishtime - starttime)
    println("Finished with ", csvfile)
end # try 

# Eerste run duurde 42.6 uur. Resultaten zijn in data/maxlik/maxlikoverallfiles.csv.gz 