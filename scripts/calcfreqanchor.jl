# calcfreq.jl
# Author: Jan van Waaij
# Created: 28 October 2022
# Last update: 20 December 2022
#
# Purpose:
# This program calculates the frequencies of the derived allele of the 
# anchor population, for each locus. 
#
# Instalation:
# Please read the install instructions in installinstructions.txt.
#
# Usage:
# open file as julia path/to/calcfreqanchor.jl
#
# NB: works only with .vcf.gz files in the folder /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/vcffiles.
#
# Output:
# output will be saved in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqanchor/filename_freqanchor.txt.gz
# The output will have a few lines of comments, each of which preceded by #, and then each line has a number between 0.0 and 1.0, which is a frequency. 

using GZip, Data

vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles")
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

vcffilespaths = joinpath.(vcffilespath, vcffiles)

vcffileswithoutextension = [file[begin:end-7] for file in vcffiles]

freqfiles = [file*"_freqanchor.txt.gz" for file in vcffileswithoutextension]
freqfilespath = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqanchor/"
freqfilespaths = joinpath.(freqfilespath, freqfiles)

@assert length(vcffilespaths) == length(freqfilespaths)
for (inputfilename, outputfilename) in zip(vcffilespaths, freqfilespaths)
    readfile = GZip.open(inputfilename)
    writefile = GZip.open(outputfilename, "w")

    # skip the header. 
    skiplines(readfile, 5);
    # save the column names. 
    columnnames = split(readline(readfile))
    numcolumns = length(columnnames)
    # The first 9 cells are information about the allele. 
    # The first exact half of the individuals are ancient individuals.
    # The second exact half of the individuals are anchor individuals. 
    lastindividualancientpopulation = 9 + (numcolumns-9)÷2 
    ancientpopulation = 10:lastindividualancientpopulation
    anchorpopulation = lastindividualancientpopulation+1:numcolumns
    # make sure that the ancientpopulation and the anchorpopulation have the same size. 
    @assert length(ancientpopulation) == length(anchorpopulation)

    comments = """
    Each line contains the frequency of the allele for the anchor population.
    original vcf file: $inputfilename
    output file: $outputfilename"""

    freq!(readfile, writefile, anchorpopulation, comments)

    close(readfile)
    close(writefile)
end

# OLD ↓

#=
filepath = ARGS[1]
filename = splitpath(filepath)[end][begin:end-7] # filename without path and file extension. 
inputfilename = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", string(filename, ".vcf.gz"))
outputfilename = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", "freqanchor", string(filename, '_', "freqanchor.txt.gz"))

using GZip, Data

readfile = GZip.open(inputfilename)
writefile = GZip.open(outputfilename, "w")

# skip the header. 
skiplines(readfile, 5);
# save the column names. 
columnnames = split(readline(readfile))
numcolumns = length(columnnames)
# The first 9 cells are information about the allele. 
# The first exact half of the individuals are ancient individuals.
# The second exact half of the individuals are anchor individuals. 
lastindividualancientpopulation = 9 + (numcolumns-9)÷2 
ancientpopulation = 10:lastindividualancientpopulation
anchorpopulation = lastindividualancientpopulation+1:numcolumns
# make sure that the ancientpopulation and the anchorpopulation have the same size. 
@assert length(ancientpopulation) == length(anchorpopulation)

comments = """
Each line contains the frequency of the allele for the anchor population.
original vcf file: $inputfilename
output file: $outputfilename"""

freq!(readfile, writefile, anchorpopulation, comments)

close(readfile)
close(writefile)
=#