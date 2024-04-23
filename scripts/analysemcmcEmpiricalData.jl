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
      CSV, # for CSV files.
 #     StatsBase, # for the mean function
#      Plots #, 
      StatsPlots,
      MCMCChains, 
      Distributed

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
@everywhere using GZip
@everywhere using DataFrames 
@everywhere using CSV
#@everywhere using StatsBase
#@everywhere using Plots 
@everywhere using StatsPlots 
@everywhere using MCMCChains 
@everywhere using Distributed

anchorpops = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, truenumind) for anchorpop in anchorpops for truenumind in truenumindvec] # 104 elementen. 4*104=416 threads

nsteps = 100_000 # number of steps in mcmc chain 
nchains = 4 # number of chains. 
trueϵ = 0.001 # true error rate used in the simulations. 

@sync @distributed for (anchorpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("EmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    plotfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "images", "empirical_data", string("EmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".png"))
    #=
    if isfile(plotfile) # plot is already made. 
        continue 
    end 
    =#
    println("Working on ", plotfile)
    df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])
    gdf = groupby(df, [:anchorpop, :trueϵ, :truenumind, :chainid])
    # chainarray = cat((Matrix(gdf[k][nsteps÷2:nsteps, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    # chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
    chainarray = cat((Matrix(gdf[k][50_000:100_000, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
    p=plot(chain, show=false);
    savefig(p, plotfile);
    #display(p)
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
    
df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])

gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :chainid])

chainarray = cat((Matrix(gdf[k][5000:6000, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
p=plot(chain)
