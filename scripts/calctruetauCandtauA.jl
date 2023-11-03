# maxlik.jl
# Author: Jan van Waaij
# Created: 2 March 2023
# Last update: 14 March 2023
#
# Purpose:
# Maximise the likelihood.
#
# Installation:
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open file as julia path/to/loglik.jl
#
# Requirement:
# Frequency files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqgenotypes
# Read files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads
# 
# Output:
# /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/likelihoods/liks3.csv
# The CSV file has 10 columns: NeAnchorPop, NeAncientPop, jointtimeyears, trueτC, trueτA, truenumind, τC, τA, numind, loglik.


using Data, # For parsing frequency and read files. 
      Likelihood, # For the likelihood function.
      GZip, # Read and write zip files. 
      DataFrames, # work with data frames.
      CSV, # convert dataframes to a CSV file. 
      Dates # to keep track of the time. 

timeatstart = now() 
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


# Testing 
# fileswithoutextensionvec = fileswithoutextensionvec[1:1]
# numberofindividualsinthesedimentvec = 1:1
# range = 1.0:1.0 

#vcffilewithoutextension = "two_demes_11000_6000_70000_0"

# sp = split(vcffilewithoutextension, '_')
# NeAncientPop = parse(Int, sp[3])
# NeAnchorPop = parse(Int, sp[4])
# jointtimeyears = parse(Int, sp[5])
# generationtime = 25
# trueτA = jointtimeyears / (2*NeAncientPop*25)
# trueτC = jointtimeyears / (2*NeAnchorPop*25)

for filewithoutextension in fileswithoutextensionvec
    # filewithoutextension is of the form two_demes_Ne1_Ne2_jointimeyears_0
    sp = split(filewithoutextension, '_')
    NeAncientPop = parse(Int, sp[3])
    NeAnchorPop = parse(Int, sp[4])
    jointtimeyears = parse(Int, sp[5])
    generationtime = 25
    trueτA = jointtimeyears / (2*NeAncientPop*generationtime)
    trueτC = jointtimeyears / (2*NeAnchorPop*generationtime)
    @show trueτC trueτA
end
