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

# trigger precompilation.
using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed

ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec] # 64 elementen. 

#randomselection = sample(1:64, 4, replace=false)

#= 
[38, 63, 36, 57]
("CHB", "YRIS", 2)
("KAR", "KARS", 4)
("CHB", "CHBS", 8)
("KAR", "CEUS", 1)

println(randomselection)
allcombinations = allcombinations[randomselection]
for combination in allcombinations
    println(combination)
end 
=# 

nodes = [
    #("node01", 1),# 20), # not used. MAX CAP FOR ME: 5 ASO MODE 6
    #("node02", 1), #12), # not used MAX CAP FOR ME: 3 ASO MODE 4
    #("node03", 24), #20), # not used MAX CAP FOR ME: 5 ASO MODE 6
    #("node04", 24), #6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    #("node05", 20), #20), #6), # 2 cores in use. MAX CAP FOR ME:5 ASO MODE 6
    #("node06", 17), # 2 core in use. MAX CAP FOR ME: 5 ASO MODE 6
    #("node07", 13), #14), # 16), # fully used MAX CAP FOR ME: 13 ASO MODE 16
    #("node08", 16), # not used MAX CAP FOR ME: 4 ASO MODE 5
    #("node11", 8), #10), # fully used MAX CAP FOR ME: 8 ASO MODE 10
    #("node12", 12), # not used MAX CAP FOR ME: 3 ASO MODE 4 
    #("node13", 12), # 2 used. MAX CAP FOR ME: 3 ASO MODE 4
    #("node14", 12), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("compute01", 60) # 50 percent used. MAX CAP FOR ME: 26 ASO MODE 32 
]

numprocs = sum(item[2] for item in nodes)

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--optimize=3", "--banner=no"] 

addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

@everywhere using Data, # For parsing frequency and read files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Distributions, # To specify the prior.
      PosteriorSampling, # To sample from the posterior.
      Dates, # for the now function 
      Base.Threads, # for threads
      Distributed


println("num procs ", nprocs())

trueϵ = 0.001 # true error rate used in the simulations. 
coverage = 64 # the coverage is 64 everywhere. 
τrange = 0.001:0.001:1.0


allcombinations2 = allcombinations[[38, 63, 36, 57]]
allcombinations3 = setdiff(allcombinations, allcombinations2)
starttime = now()
@sync @distributed for (anchorpop, ancientpop, truenumind) in allcombinations3    
    starttime = now()
    frequencyfile = "data/frequencies/freqanchor/Madsmodel_freq_anchor_$(anchorpop).txt.gz"
    if !isfile(frequencyfile)
        println(frequencyfile, " does not exist.")
        continue
    end 
    freqs = parsefreqfile(frequencyfile)
    readsfile = string("data/reads/Mads_model_pop_", ancientpop, "_num_ind_",
    truenumind, "_cov_", coverage, "_reads_with_errors.txt.gz")
    if !isfile(readsfile)
        println(readsfile, " does not exist.")
        continue
    end 
    derivedreadsvec, coveragevec = parsereadfile(readsfile)
    @assert length(unique(coveragevec)) == 1
    @assert coveragevec[1] == coverage
    @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "exactposterior", string("Madsmodel_ancient_pop_exact_posterior_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    #=
    if isfile(zippedcsvfile)
        println(zippedcsvfile, " does already exist. Do not overwrite existing files!")
        continue 
    end
    =# #pid 614428
    df = DataFrame(anchorpop = String[],
                   ancientpop = String[],
                   trueϵ = Float64[],
                   truenumind = Int[],
                   coverage = Int[],
                   τCvalue = Float64[],
                   τAvalue = Float64[],
                   logjointprob = Float64[]
                  )
    CSV.write(csvfile, df)
    results = ExactPosterior(truenumind, τrange, trueϵ, coverage, derivedreadsvec, freqs, 1000)
    ncalculations = size(results, 2)
    local df = DataFrame(anchorpop = repeat([anchorpop], ncalculations),
                         ancientpop = repeat([ancientpop], ncalculations),
                         trueϵ = repeat([trueϵ], ncalculations),
                         truenumind = repeat([truenumind], ncalculations),
                         coverage = repeat([coverage], ncalculations),
                         τCvalue = results[1,:],
                         τAvalue = results[2,:],
                         logjointprob = results[3,:]
                        )
    CSV.write(csvfile, df, append=true)
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
finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("TIME TO FINISH everything ", numhours, " and ", remainingtimeminutes, '.')

#=
Finished duurde ongeveer 30-38 uur. 
=#