# analyseposterior.jl
# Author: Jan van Waaij
# Created: 9 May 2024
# Last update: 9 May 2024
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

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      CSV, # for CSV files.
      StatsBase, # for the mean function
      Plots, # plotting. 
      StatsPlots, 
      MCMCChains, # plotting of mcmc chains. 
      Distributed, # distributed programming
      Latexify # make latex code from Julia objects. 

numberofindividualsinthesedimentvec = 1:8
coveragesvec = 1:100

nsteps = 100_000 # number of steps in mcmc chain 
nchains = 4 # number of chains. 
trueϵ = 0.003 # true error rate used in the simulations. 

results = DataFrame(
    truenumind = Int[],
    coverage = Int[], 
    n_low = Float64[],
    n_median = Float64[],
    n_high = Float64[],
    τC_low = Float64[],
    τC_median = Float64[],
    τC_high = Float64[],
    τA_low = Float64[],
    τA_median = Float64[],
    τA_high = Float64[],
    ϵ_low = Float64[],
    ϵ_median = Float64[],
    ϵ_high = Float64[],
)

allind = 1:8
allcov = 1:100
allcombinations = [(numind, cov) for numind in allind for cov in allcov]

for (numberofindividualsinthesediment, coverage) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", coverage, "_MCMC.csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame, header = ["trueϵ", "truenumind", "coverage", "nsample", "τCsample", "τAsample", "ϵsample", "accepted", "logjointprob", "chainid"])
    if eltype(df.trueϵ) <: AbstractString
        df = CSV.read(zippedcsvfile, DataFrame)
    end 
    df = df[end÷2:10:end,:]
    gdf = groupby(df, [:truenumind, :coverage, :trueϵ, :chainid])
    n_low, n_median, n_high = quantile(df.nsample, [0.025, 0.5, 0.975])
    τC_low, τC_median, τC_high = quantile(df.τCsample, [0.025, 0.5, 0.975])
    τA_low, τA_median, τA_high = quantile(df.τAsample, [0.025, 0.5, 0.975])
    ϵ_low, ϵ_median, ϵ_high = quantile(df.ϵsample, [0.025, 0.5, 0.975])
    push!(results, (truenumind = numberofindividualsinthesediment, 
                    coverage = coverage, 
                    n_low = n_low, 
                    n_median = n_median, 
                    n_high = n_high, 
                    τC_low = τC_low, 
                    τC_median = τC_median, 
                    τC_high = τC_high, 
                    τA_low = τA_low, 
                    τA_median = τA_median, 
                    τA_high = τA_high, 
                    ϵ_low = ϵ_low, 
                    ϵ_median = ϵ_median, 
                    ϵ_high = ϵ_high))
end 

results_sorted = sort(results, [:truenumind, :coverage])

function afronden(x::Float64) 
    if isinteger(x)
        return string(Int(x))
    else
        y = round(x, sigdigits=3)
        if isinteger(y)
            return string(Int(y))
        else
            return y 
        end
    end 
end 
afronden(x::Int) = string(x)
afronden(x) = x 

results_afgerond_sorted = afronden.(results_sorted)

for numind in numberofindividualsinthesedimentvec
    results_sorted_fixed_n = results_sorted[results_sorted.truenumind .== numind, :]
    local coveragesvec = sort(unique(results_sorted_fixed_n.coverage), rev=true)
    singletonconfidenceinterval = Int[]
    highestcoverageforwhichintervalisnotasingleton = 0 
    for coverage in coveragesvec 
        results_sorted_fixed_n_and_fixed_coverage = results_sorted_fixed_n[results_sorted_fixed_n.coverage .== coverage, :]
        @assert size(results_sorted_fixed_n_and_fixed_coverage, 1) == 1 
        if results_sorted_fixed_n_and_fixed_coverage[1, :n_low] == results_sorted_fixed_n_and_fixed_coverage[1, :n_high]
            push!(singletonconfidenceinterval, coverage)
        else
            highestcoverageforwhichintervalisnotasingleton = coverage
            break # break the inner loop. 
        end 
    end 
    latexcode = latexify(results_afgerond_sorted[(results_afgerond_sorted.truenumind .== string(numind)) .& (results_sorted.coverage .≤ highestcoverageforwhichintervalisnotasingleton), [:truenumind, :coverage, :n_low, :n_median, :n_high]], env=:table)
    filename = "data/latexsnippets/confidence_interval_n_num_ind_$(numind).tex"
    openedfile = open(filename, "w")
    write(openedfile, "\\begin{table}[H]\n")
    write(openedfile, latexcode)
    write(openedfile, "\\end{table}\n ")
    close(openedfile)
    
    singletoncrediblesetphrase = "For coverages "
    sort!(singletonconfidenceinterval)
    for coverage in singletonconfidenceinterval[begin:end-1]
        singletoncrediblesetphrase *= string(coverage) * ", "
    end 
    singletoncrediblesetphrase = singletoncrediblesetphrase * "and $(singletonconfidenceinterval[end]) the credible sets are \\{$(numind)\\}. "
    filename = "data/latexsnippets/singleton_interval_phrase_num_ind_$(numind).tex"
    write(filename, singletoncrediblesetphrase)
    
    latexcode = latexify(results_afgerond_sorted[results_afgerond_sorted.truenumind .== string(numind), [:truenumind, :coverage, :τC_low, :τC_median, :τC_high]], env=:table)
    filename = "data/latexsnippets/confidence_interval_tauC_num_ind_$(numind).tex"
    openedfile = open(filename, "w")
    write(openedfile, "\\begin{table}[H]\n")
    write(openedfile, latexcode)
    write(openedfile, "\\end{table}\n ")
    close(openedfile)
    
    latexcode = latexify(results_afgerond_sorted[results_afgerond_sorted.truenumind .== string(numind), [:truenumind, :coverage, :τA_low, :τA_median, :τA_high, :ϵ_low, :ϵ_median, :ϵ_high ]], env=:table)
    filename = "data/latexsnippets/confidence_interval_tauA_and_eps_num_ind_$(numind).tex"
    openedfile = open(filename, "w")
    write(openedfile, "\\begin{table}[H]\n")
    write(openedfile, latexcode)
    write(openedfile, "\\end{table}\n ")
    close(openedfile)
end 

