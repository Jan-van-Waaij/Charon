# analyseposterior.jl
# Author: Jan van Waaij
# Created: 6 March 2024
# Last update: 6 March 2024
#
# Purpose:
# This is a collection of Julia commands to analyse the MCMC samples.
#
# Installation:
# You need the GZip, DataFrames, Plots and CSV packages. 
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open a Julia terminal and use snippets of the code interactively. 
#
# Requirement:
# a csv.gz file with MCMC samples from /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/mcmcsamples

# UNKNOWN N and with errors:

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      CSV, 
      CairoMakie,
      StatsBase 

folder = "/home/projects/gabriel/Jan_allelefreq/hg38/bams/"
files = readdir(folder)
freqfiles = files[contains.(files, r"IBS_ind\d_freq_[A-Z]{3}.gz")]
populations = unique([freqfile[end-5:end-3] for freqfile in freqfiles])
allthecombinations = [(population, numind) for population in populations for numind in [1, 2, 4, 8]]
    
dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/images/empirical_data/"

nsteps = 100_000 # number of steps in mcmc chain 
nchains = 4 # number of chains. 

results = DataFrame(
    numind = Int[],
    ancient = String[],
    anchor = String[], 
    n_low = Float64[],
    n_median = Float64[],
    n_high = Float64[],
    τC_low = Float64[],
    τC_median = Float64[],
    τC_high = Float64[],
    τA_low = Float64[],
    τA_median = Float64[],
    τA_high = Float64[],
    τCplusτA_low = Float64[], 
    τCplusτA_median = Float64[], 
    τCplusτA_high = Float64[],
    maxτCτA_low = Float64[], 
    maxτCτA_median = Float64[], 
    maxτCτA_high = Float64[],
    ϵ_low = Float64[],
    ϵ_median = Float64[],
    ϵ_high = Float64[],
)
  
csvfilefolder = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "EmpiricalData")

for (population, numind) in allthecombinations
    csvfile = joinpath(csvfilefolder, string("EmpiricalData", "_anchor_pop_", population, "_true_num_ind_", numind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        printstyled(string(zippedcsvfile, " is missing.\n"), color = :red)
        continue 
    end 
    
    df = CSV.read(zippedcsvfile, DataFrame)
    df = df[end÷2:10:end,:]
    n_low, n_median, n_high = quantile(df.nsample, [0.025, 0.5, 0.975])
    τC_low, τC_median, τC_high = quantile(df.τCsample, [0.025, 0.5, 0.975])
    τA_low, τA_median, τA_high = quantile(df.τAsample, [0.025, 0.5, 0.975])
    ϵ_low, ϵ_median, ϵ_high = quantile(df.ϵsample, [0.025, 0.5, 0.975])
    τCplusτA_low, τCplusτA_median, τCplusτA_high = quantile(df.τCsample+df.τAsample, [0.025, 0.5, 0.975])
    maxτCτA_low, maxτCτA_median, maxτCτA_high = quantile(max.(df.τCsample, df.τAsample), [0.025, 0.5, 0.975])
    push!(results, (numind = numind,
                    ancient = "IBS",
                    anchor = population, 
                    n_low = n_low,
                    n_median = n_median,
                    n_high = n_high,
                    τC_low = τC_low,
                    τC_median = τC_median,
                    τC_high = τC_high,
                    τA_low = τA_low,
                    τA_median = τA_median,
                    τA_high = τA_high,
                    τCplusτA_low = τCplusτA_low, 
                    τCplusτA_median = τCplusτA_median, 
                    τCplusτA_high = τCplusτA_high,
                    maxτCτA_low = maxτCτA_low, 
                    maxτCτA_median = maxτCτA_median, 
                    maxτCτA_high = maxτCτA_high,
                    ϵ_low = ϵ_low,
                    ϵ_median = ϵ_median,
                    ϵ_high = ϵ_high,))
end 

filename = "resultsEmpirical.csv" 
CSV.write(filename, results)
cmd = `gzip -f $filename`
run(cmd)

matrices = Vector{Matrix{String3}}(undef, 4)
numindvec = [1, 2, 4, 8]
sortedresults = DataFrame(
method = String[],
numind = Int[],
ord1 = String[],
ord2 = String[],
ord3 = String[],
ord4 = String[],
ord5 = String[],
ord6 = String[],
ord7 = String[],
ord8 = String[],
ord9 = String[],
ord10 = String[],
ord11 = String[],
ord12 = String[],
ord13 = String[],
ord14 = String[],
ord15 = String[],
ord16 = String[],
ord17 = String[],
ord18 = String[],
ord19 = String[],
ord20 = String[],
ord21 = String[],
ord22 = String[],
ord23 = String[],
ord24 = String[],
ord25 = String[],
ord26 = String[],
)

for index in eachindex(numindvec)
    k = numindvec[index]
    resultsf = filter(:numind => x -> x == k, results)
    τCest = resultsf.:τC_median
    τAest = resultsf.:τA_median

    permτC = sortperm(τCest)

    permτA = sortperm(τAest)

    permτCplusτA = sortperm(τCest+τAest)

    permmaxτCτA = sortperm(max.(τCest, τAest))

    pops = resultsf.:anchor

    sortedpopsbyτC = pops[permτC]
    sortedpopsbyτA = pops[permτA]
    sortedpopsbyτCplusτA = pops[permτCplusτA]
    sortedpopsbymaxτCτA = pops[permmaxτCτA]

    push!(sortedresults, ["τC"; k; sortedpopsbyτC])
    push!(sortedresults, ["τA"; k; sortedpopsbyτA])
    push!(sortedresults, ["τCplusτA"; k; sortedpopsbyτCplusτA])
    push!(sortedresults, ["maxτCτA"; k; sortedpopsbymaxτCτA])
end 

k = 2
resultsf = filter(:numind => x -> x == k, results)

τCest = resultsf.:τC_median
permτC = sortperm(τCest)

pops = resultsf.:anchor

ys = sortedpopsbyτC = pops[permτC]

lows = resultsf.:τC_low[permτC]
xs = medians = resultsf.:τC_median[permτC]
highs = resultsf.:τC_high[permτC]

lowerrors = medians - lows 
higherrors = highs - medians

f = Figure(; size=(1600, 1600))

ax = Axis(f[1,1], yticks = (1:26, ys))

errorbars!(ax, xs, 1:26, lowerrors, higherrors; whiskerwidth = 10, color = :magenta2, direction = :x)


methods = ["τC", "τA", "τCplusτA", "maxτCτA"]
labels = [L"\tau_C", L"\tau_A", L"\tau_C+\tau_A", L"\max(\tau_C\,,\, \tau_A)"]
filenameparts = ["tauC", "tauA", "tauCplustauA", "maxtauCtauA"]
numindvec = [1, 2, 4, 8]

f = Figure(; size=(2400, 2400))

colours = [
    :red :turquoise1 :maroon1 :chocolate 
    :magenta2 :darkgreen :grey63 :navajowhite1  
    :cyan :yellow3 :darkorange1 :red4 
    :lightsalmon :yellow1 :green3 :teal 
]

axes = [Axis(f[numindindex, methodindex], title = "Number of individuals $(numindvec[numindindex])", xlabel = labels[methodindex], xlabelsize=30) for numindindex in 1:4, methodindex in 1:4]

linkxaxes!(axes[:, 1]...)
linkxaxes!(axes[:, 2]...)
linkxaxes!(axes[:, 3]...)
linkxaxes!(axes[:, 4]...)


for methodindex in eachindex(methods, labels)
    method = methods[methodindex] 
    label = labels[methodindex]
    for numindindex in eachindex(numindvec)
        numind = numindvec[numindindex]
        resultsf = filter(:numind => x -> x == numind, results)
        if method == "τC"
            est = resultsf.:τC_median
            lows = resultsf.:τC_low
            highs = resultsf.:τC_high
        elseif method == "τA"
            est = resultsf.:τA_median
            lows = resultsf.:τA_low
            highs = resultsf.:τA_high
        elseif method == "τCplusτA"
            est = resultsf.:τCplusτA_median
            lows = resultsf.:τCplusτA_low
            highs = resultsf.:τCplusτA_high
        elseif method == "maxτCτA"
            est = resultsf.:maxτCτA_median
            lows = resultsf.:maxτCτA_low
            highs = resultsf.:maxτCτA_high
        else
            error("something went wrong.")
        end 
        
        lowerrors = est - lows 
        higherrors = highs - est 
        
        perm = sortperm(est)
        
        xs = est[perm]
        lowerrors = lowerrors[perm]
        higherrors = higherrors[perm]

        pops = resultsf.:anchor[perm]
        ys = 1:length(pops)

        ax = axes[numindindex, methodindex] #Axis(f[numindindex, methodindex], yticks = (ys, pops), title = "Number of individuals $numind", xlabel = label, xlabelsize=30)
        ax.yticks = (ys, pops)
        errorbars!(ax, xs, ys, lowerrors, higherrors; whiskerwidth = 10, color = colours[numindindex, methodindex], direction = :x)
        scatter!(ax, xs, ys; markersize = 5, color = :black)
    end 
end 

save("data/images/empirical_data/credibleintervalsempircaldatasortedeverything.png", f)

for methodindex in eachindex(methods, labels, filenameparts)
    method = methods[methodindex] 
    label = labels[methodindex]
    for numindindex in eachindex(numindvec)
        numind = numindvec[numindindex]
        resultsf = filter(:numind => x -> x == numind, results)
        if method == "τC"
            est = resultsf.:τC_median
            lows = resultsf.:τC_low
            highs = resultsf.:τC_high
        elseif method == "τA"
            est = resultsf.:τA_median
            lows = resultsf.:τA_low
            highs = resultsf.:τA_high
        elseif method == "τCplusτA"
            est = resultsf.:τCplusτA_median
            lows = resultsf.:τCplusτA_low
            highs = resultsf.:τCplusτA_high
        elseif method == "maxτCτA"
            est = resultsf.:maxτCτA_median
            lows = resultsf.:maxτCτA_low
            highs = resultsf.:maxτCτA_high
        else
            error("something went wrong.")
        end 
        
        lowerrors = est - lows 
        higherrors = highs - est 
        
        perm = sortperm(est)
        
        xs = est[perm]
        lowerrors = lowerrors[perm]
        higherrors = higherrors[perm]

        pops = resultsf.:anchor[perm]
        ys = 1:length(pops)

        f = Figure(; size=(1200, 1200))
        ax = Axis(f[1,1], yticks = (ys, pops), title = "Number of individuals $numind", xlabel = label, xlabelsize=30)
        errorbars!(ax, xs, ys, lowerrors, higherrors; whiskerwidth = 10, color = colours[numindindex, methodindex], direction = :x)
        scatter!(ax, xs, ys; markersize = 5, color = :black)
        save("data/images/empirical_data/credibleintervalsempircaldatasorted$(filenameparts[methodindex])numind$(numind).png", f)
    end 
end 


f = Figure(; size=(3200, 800))

axes = [Axis(f[1, a], xlabel = "number of intervals", title = "Credible intervals for the number of intervals") for a in 1:4]

numindvec = [1, 2, 4, 8]
for (numindindex, numind) in pairs(numindvec)
    resultsf = filter(:numind => t -> t == numind, results)
    #sort!(resultsf, :τCplusτA_median)

    xs = resultsf.:n_median
    higherrors = resultsf.:n_high - xs 
    lowerrors = xs - resultsf.:n_low

    pops = resultsf.:anchor 

    ys = 1:size(resultsf, 1)

    yticks = (ys, pops)
    
    ax = axes[numindindex]
    #ax.yticks = yticks  

    vlines!(ax, [numind], color=:blue , linewidth = 0.5)

    errorbars!(ax, xs, ys, lowerrors, higherrors; whiskerwidth = 10, color = :red, direction = :x)

end 

f 



f = Figure(; size=(800, 800))

ax = Axis(f[1, 1], xlabel = "number of intervals", title = "Credible intervals for the number of intervals")

numind = 4
resultsf = filter(:numind => t -> t == numind, results)
    #sort!(resultsf, :τCplusτA_median)

xs = resultsf.:n_median
higherrors = resultsf.:n_high - xs 
lowerrors = xs - resultsf.:n_low

pops = resultsf.:anchor 

ys = 1:size(resultsf, 1)

yticks = (ys, pops)
    
ax.yticks = yticks  

#errorbars!(ax, xs, ys, lowerrors, higherrors; whiskerwidth = 10, color = :red, direction = :x)
vlines!(ax, [4.5], ymin=minimum(ys), ymax = maximum(ys), color=:blue , linewidth = 20)


f


f = Figure(; size = (800, 800))
ax = Axis(f[1,1])
vlines!(ax, 4.0)
