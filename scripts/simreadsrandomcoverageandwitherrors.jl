# simreads.jl
# Author: Jan van Waaij
# Created: 19 September 2023
# Last update: 19 September 2023
#
# Purpose:
# Simulate reads for the numberofindividualsinthesediment individuals in the sediment, 
# who are randomly chosen from the ancient population in filename.vcf.gz. The coverage is set with coverage. 
#
# Installation:
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open file as julia path/to/simreads.jl path/to/filename.vcf.gz numberofindividualsinthesediment coverage
#
# NB: works only with .vcf.gz files in the folder /home/projects/SedaDNA_Demo/simulationsMSPRIME/data.
#
# Output:
# output will be saved in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads/filename_numberofindividualsinthesediment_coverage_reads.txt.gz
# The output will have a few lines of comments, each of which preceded by #, and then each line has two integers, of which the first is the number
# of derived reads and the second is the total number of reads (a.k.a the coverage).

using Data, GZip, StatsBase, Distributions

vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles")
readfilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads")
allfilesunfiltered = readdir(vcffilespath)
vcffiles = String[]
for file in allfilesunfiltered 
    if file[end-6:end] != ".vcf.gz" # remove non VCF files.
        continue
    end 
    sp = split(file, '_')
    if sp[3] == "1000" # remove two_demes_1000_10_0.vcf.gz
        continue 
    end 
    push!(vcffiles, file)
end 

vcffileswithoutextensionvec = [file[begin:end-7] for file in vcffiles]

numberofindividualsinthesedimentvec = 1:8 # [1, 2, 4, 8]
ϵ = 0.003

coverage = Poisson(1)

for vcffilewithoutextension in vcffileswithoutextensionvec
    # Here I determine which individuals are in the ancient population
    # and which are in the anchor population. 
    vcffilepath = joinpath(vcffilespath, vcffilewithoutextension * ".vcf.gz")
    openedvcffile = GZip.open(vcffilepath)
    # skip the header. 
    skiplines(openedvcffile, 5)
    # save the column names. 
    columnnames = split(readline(openedvcffile))
    close(openedvcffile)
    numcolumns = length(columnnames)
    lastindividualancientpopulation = 9 + (numcolumns-9)÷2 
    ancientpopulation = 10:lastindividualancientpopulation
    anchorpopulation = lastindividualancientpopulation+1:numcolumns
    @assert length(ancientpopulation) == length(anchorpopulation)    
    
    for numberofindividualsinthesediment in numberofindividualsinthesedimentvec
        # I choose numberofindividualsinthesediment individuals from the ancient population. 
        individualsinthesediment = sample(ancientpopulation, numberofindividualsinthesediment, 
                                            replace=false, ordered=true)
        # I open the file where I while store the reads. 
        readfile = string(vcffilewithoutextension, '_', 
            numberofindividualsinthesediment, '_', "1X_with_errors.txt.gz")
        pathtoreadfile = joinpath(readfilespath, readfile)
        openedreadfile = GZip.open(pathtoreadfile, "w")

        openedvcffile = GZip.open(vcffilepath)
        skiplines(openedvcffile, 6)
        
        comments = """
        Each line contains the number of derived reads and the total number of reads, seperated by a tab.  
        original vcf file: $vcffilepath 
        Individuals in the sediment: $individualsinthesediment
        Total number of individuals in the sediment: $numberofindividualsinthesediment
        coverage: $coverage"""
        reads(openedvcffile, openedreadfile, individualsinthesediment, coverage, comments, ϵ)
        close(openedvcffile)
        close(openedreadfile)
    end  
    close(openedvcffile)
end






























#=

filepath = ARGS[1] 
numberofindividualsinthesediment = parse(Int, ARGS[2])
coverage = parse(Int, ARGS[3])

filepath = ARGS[1]
filename = splitpath(filepath)[end][begin:end-7] # filename without path and file extension. 
inputfilename = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles", string(filename, ".vcf.gz"))
outputfilename = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", string(filename, '_'
    , numberofindividualsinthesediment, '_', coverage, '_' ,"reads.txt.gz"))

readfile = GZip.open(inputfilename)
writefile = GZip.open(outputfilename, "w")

# skip the header. 
skiplines(readfile, 5);
# save the column names. 
columnnames = split(readline(readfile))
numcolumns = length(columnnames)
lastindividualancientpopulation = 9 + (numcolumns-9)÷2 
ancientpopulation = 10:lastindividualancientpopulation
anchorpopulation = lastindividualancientpopulation+1:numcolumns
@assert length(ancientpopulation) == length(anchorpopulation)

# select numberofindividualsinthesediment individuals in population1.
individualsinthesediment = sample(ancientpopulation, numberofindividualsinthesediment, replace=false, ordered=true)

comments = """
Each line contains the number of derived reads and the total number of reads, seperated by a tab.  
original vcf file: $inputfilename 
Individuals in the sediment: $individualsinthesediment
Total number of individuals in the sediment: $numberofindividualsinthesediment
coverage: $coverage"""

reads(readfile, writefile, individualsinthesediment, coverage, comments)

close(readfile)
close(writefile)

=#