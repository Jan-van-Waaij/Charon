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

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/images/empirical_data/"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--optimize=3"]
      
nodes = [
   #= ("compute01", 30), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32 
    ("node05", 6), #6), # 2 cores in use. MAX CAP FOR ME:5 ASO MODE 6
    ("node03", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node04", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node12", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4 
    ("node01", 6), # not used. MAX CAP FOR ME: 5 ASO MODE 6
    ("node02", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("node06", 6), # 2 core in use. MAX CAP FOR ME: 5 ASO MODE 6
    ("node07", 13), #14), # 16), # fully used MAX CAP FOR ME: 13 ASO MODE 16
    ("node13", 4), # 2 used. MAX CAP FOR ME: 3 ASO MODE 4
    ("node11", 10), #10), # fully used MAX CAP FOR ME: 8 ASO MODE 10
  =#
    ("node14", 16) #, # not used MAX CAP FOR ME: 3 ASO MODE 4
#    ("node08", 5), # not used MAX CAP FOR ME: 4 ASO MODE 5
]
sum(node[2] for node in nodes)


addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

@everywhere ENV["GKSwstype"]="100"
@everywhere using GZip, # to work with zipped files,
DataFrames, # to work with DataFrames,
CSV, # for CSV files.
StatsBase, # for the mean function
Plots, # plotting. 
StatsPlots, 
MCMCChains, # plotting of mcmc chains. 
Distributed # distributed programming

number_of_individuals_vec = 1:8
coverage_vec = 2:20 
allcombinations = [(numind, coverage) for numind in number_of_individuals_vec for coverage in coverage_vec]; # 152 elements. 

nsteps = 100_000 # number of steps in mcmc chain 
nchains = 4 # number of chains. 
trueϵ = 0.003 # true error rate used in the simulations. 

@sync @distributed for (numind, coverage) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("knowntau=0.25_num_ind_", numind, "_cov_", coverage, "_standard_MCMC.csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    plotfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "images", "knowntau", string("knowntau=0.25_num_ind_", numind, "_cov_", coverage, "_standard_MCMC.png"))
    if isfile(plotfile) # plot is already made. 
        continue 
    end 
    println("Working on ", plotfile)
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:truenumind, :coverage, :trueϵ, :chainid])
    chainarray = cat((Matrix(gdf[k][50_000:100_000, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
    p=plot(chain, show=false);
    savefig(p, plotfile);
end 

tableτA = Array{Vector{Float64}}(undef, 19, 8)
tableτC = Array{Vector{Float64}}(undef, 19, 8)
# table = zeros(19,8)

for (numind, coverage) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("knowntau=0.25_num_ind_", numind, "_cov_", coverage, "_standard_MCMC.csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:truenumind, :coverage, :trueϵ, :chainid])
    df = vcat((gdf[k][end÷2+1:10:end, :] for k in eachindex(gdf))...)
    tableτC[coverage-1, numind] = round.(quantile(df.:τCsample, [0.025, 0.5, 0.975]), digits=2)
    tableτA[coverage-1, numind] = round.(quantile(df.:τAsample, [0.025, 0.5, 0.975]), digits=2)
    #table[coverage-1, numind] = mean(df.τAsample)
end 

table = round.(table, digits=3)
latexify(table)

tabletuplesτC = similar(tableτC, Tuple{Float64, Float64, Float64})

for index in eachindex(tableτC, tabletuples)
    tabletuplesτC[index] = tuple(tableτC[index]...)
end 

tableτCstrings = string.(tableτC)
latexify(tableτCstrings)

tabletuplesτA = similar(tableτA, Tuple{Float64, Float64, Float64})

for index in eachindex(tableτA, tabletuplesτA)
    tabletuplesτA[index] = tuple(tableτA[index]...)
end 

tableτAstrings = string.(tabletuplesτA)
s = latexify(tableτAstrings[:, 5:8])

sp = split(s, '\n')
coverage = 2 
news = ""
for index in eachindex(sp)[4:end-4]
    news = news * string(coverage) * " & " * sp[index] * "\n\\hline\n"
    global coverage += 1 
end 
print(news)




tablecontained = similar(table, Bool)
for index in eachindex(table, tablecontained)
    tablecontained[index] = table[index][1] ≤ 0.25 ≤ table[index][3]
end



function printresults(anchorpop, ancientpop, truenumind, chainid, rate)
    output = string("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " chain ", chainid)
    output = output * string('\n', rate, "\n")
    println(output)
end 
# calculate acceptancerate
for (anchorpop, ancientpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("simprop_Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    for df in gdf 
        rate = mean(df.:accepted[50_000:100_000])
        printresults(anchorpop, ancientpop, truenumind, df.:chainid[1], rate)
    end 
end 

anchorpop = "ACB"
truenumind = 2

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("EmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
zippedcsvfile = csvfile * ".gz"
    
zippedcsvfile = "data/mcmcsamples/storeliksEmpericalData_anchor_pop_MSL_true_num_ind_8.csv"
df = CSV.read(zippedcsvfile, DataFrame)

df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])

gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :chainid])

chainarray = cat((Matrix(gdf[k][5000:6000, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
p=plot(chain)
