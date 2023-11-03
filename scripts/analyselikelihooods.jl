# analyselikelihoods.jl
# Author: Jan van Waaij
# Created: 14 March 2023
# Last update: 
#
# Purpose:
# This is a collection of Julia commands to analyse the likelihoods, especially
# determining the maximum likelihood estimators. 
#
# Installation:
# You need the GZip, DataFrames, Plots and CSV packages. 
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open a Julia terminal and use snippets of the code interactively. 
#
# Requirement:
# a csv.gz file from /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/likelihoods

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      Plots, # for plotting,
      CSV # for CSV files. 

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "mcmcnolocalmaximacoverage100.csv.gz")

# openedcsvfile = GZip.open(csvfile)
df = CSV.read(csvfile, DataFrame)
# close(openedcsvfile)

gdf = groupby(df, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :truenumind])


cdf = combine(gdf, [:loglik, :numind, :τC, :τA] 
                    => ((x,y...) -> reshape([z[argmax(x)] for z in y], 1, :)) 
                    => [:maxliknumind, :maxlikτC, :maxlikτA])

# shuffle the columns. 
scdf = cdf[:, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :truenumind, :maxliknumind, :trueτC, :maxlikτC, :trueτA, :maxlikτA]]

# show all columns. 
show(scdf, allrows=true, eltypes=false)