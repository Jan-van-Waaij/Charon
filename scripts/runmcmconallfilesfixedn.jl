# runmcmconallfiles.jl
# Author: Jan van Waaij
# Created: 4 May 2023
# Last update: 12 June 2023
#
# Purpose:
# Run MCMC.
# 
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia -O3 --threads=auto path/to/runmcmconallfiles.jl
# -O3 finds the most optimal optimisation. 
# --threads=auto, let Julia figure out the number of threads. 
# --threads=5 uses 5 threads. 
# If threads is not given, Julia uses one thread, and the chains cannot be 
# simulated simultaniously. 
# 
# Requirement:
# Frequency files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqgenotypes
# Read files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads
# 
# Output:
# Zipped csv File in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/mcmcsamples
# The CSV file has 10 columns: NeAnchorPop, NeAncientPop, jointtimeyears, trueτC, trueτA, truenumind, τC, τA, numind, loglik, chainid.


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

numberofindividualsinthesedimentvec = [4, 8]#, 16, 32, 64] # [1, 2, 4]#, 8, 16, 32, 64]
#coveragevec = [20, 100]
maxnumreads = [5, 10, 15, 20]#, 100]
# 16 individuals and coverage 20. gaat het mis.
# kijk of het zonder threads wel goed gaat.

nchains = 4 # number of chains


csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "allsamplesfixedndifferentcoveragesanddifferentnumindividualsmanyindividuals.csv.gz")
if isfile(csvfile)
    println(csvfile, " does already exist. Do not overwrite existing files!")
    exit(1)
end
#openedcsvfile = GZip.open(csvfile, "a") # the a of append. 

df = DataFrame(
                    NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    trueτC = Float64[],
                    trueτA = Float64[],
                    truenumind = Int[],
                    coverage = Int[],
                    τCsample = Float64[],
                    τAsample = Float64[],
                    accepted = Bool[],
                    chainid = Int[]
               )

CSV.write(csvfile, df)

nsteps = 100_000 
prior = Product([Uniform(), Uniform()]) # uiform prior on the interval 
# [0,1] on both τC and τA.

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
        NeAncientPop = parse(Int, sp[3]) # convert strings to integers. 
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
            for maxreads in maxnumreads
                @show filewithoutextension truenumind maxreads
                readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", 
                    string(filewithoutextension, '_', truenumind, '_', maxreads, "_reads.txt.gz") )
                if !isfile(readsfile)
                    println(readsfile, " does not exist.")
                    continue
                end 
                gzopen(readsfile) do f 
                    global derivedreadsvec, coveragevec = parsereadfile(f)
                end
                @assert length(unique(coveragevec)) == 1 # make sure 
                # coveragevec is a constant vector.
                coverage = coveragevec[1]
                @assert maxreads == coverage
                @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
                
                #init = (τC = rand(), τA = rand())

                results = MHsampler(nchains, nsteps, prior, coverage, derivedreadsvec, freqs, truenumind, 25_000) 
                # results is a vector with triples containing vectors.

                for i in eachindex(results) 
                    τCsample, τAsample, acceptedvec = results[i]
                    chainid = repeat([Int(i)], nsteps)
                    local df = DataFrame(NeAnchorPop = repeat([NeAnchorPop], nsteps),
                                NeAncientPop = repeat([NeAncientPop], nsteps),
                                jointtimeyears = repeat([jointtimeyears], nsteps),
                                trueτC = repeat([trueτC], nsteps),
                                trueτA = repeat([trueτA], nsteps),
                                truenumind = repeat([truenumind], nsteps),
                                coverage = repeat([maxreads], nsteps),
                                τCsample = τCsample,
                                τAsample = τAsample,
                                accepted = acceptedvec,
                                chainid = chainid)
                    CSV.write(csvfile, df, append=true)
                end 

                println("finished ", filewithoutextension, " true num ind ", truenumind)
            end 
        end 
    end # for filewithoutextension
finally
    #close(openedcsvfile)
    finishtime = now()
    println(finishtime - starttime)
    println("everything written to ")
    println("finished with ", csvfile)
end # try 

# old version 96459094 miliseconds, approximately 26.8 uur. File is allsamples4. new version 4.4 uur. Dat is 6.1 keer sneller dat is allsamples5. 
# faster version for uniform: 4.36 uur, dus nauwelijks beter.
# Inline version: 4.50 uur. 
# Conclusie: laat alles zoals het is.
# Met Eigen decomposition: 4.20 uur. Ietjes beter, maar
# niet zo goed als dat ik gehoopt had.
# Met @inline: 4.22 uur maakt dus geen kloot uit. 
# Met -O3 optie: 4.38 uur. 
# Met StaticArrays en optie -O3: 2.61 uur WOW Dus zo'n 60 sneller.
# Met StaticArrays zonder optie: 2.62 uur maakt dus niet uit.
# Met StaticArrays, -O3 en checkbounds=no. 2.61 uur, maakt dus ook niets uit. 
# Het is dus StaticArrays dat de magic doet.
# Alles met static arrays: dus ook hvals, binomialprobs, en d. 
# Ik denk dat hvals en d zo groot zijn, dat dit misschien wel te gek is. We moeten 
# kijken. Eventueel deze twee met SizedArray.Stackoverflow error =) 
# Met SizedArrays voor grote arrays: 1.93 uur. Dus weer 35 % sneller.
# Same but now d normal array: 1.96 uur.
# Same but now d StaticArray: horrible still running after two days. Killed it. 
# binomialprobs and d sizedarrays: 2.08 dus slechter. 
# Optimale setting 1.93 uur. 
# Nu Likelihood and PosteriorSampling opschonen.[]
# new mcmcsampler. Nieuwe proposal van de vorm Normal(Σ) en zonder uniforme prior aan te nemen. 1.89 uur. Met multitreading: 2.05.

# Ik wil testen wat de invloed van het aantal individuen is. In onze vorige setup toetsten we 1, 2, en 4 (3, 5, en 6 bestaan niet.)
# Nu toets is met [4, 4, 4]. Als het veel uitmaakt, moet dit significant veel langzamer zijn. nu 2.91 uur. Met multitreading: 3.00 uur. 

