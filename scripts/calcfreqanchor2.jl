# calcfreq.jl
# Author: Jan van Waaij
# Created: 2 November 2023
# Last update: 2 November 2023
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

inputfilename = "data/vcffiles_demes/two_demes_with_demes_tauC=0.25tauA=0.25numgen=1000effpop=200012.vcf.gz"
outputfilename = "data/frequencies/freqanchor/two_demes_demes_freqanchor.txt.gz"

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

freq(readfile, writefile, anchorpopulation, comments)

close(readfile)
close(writefile)
