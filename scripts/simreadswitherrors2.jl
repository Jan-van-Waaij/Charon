# simreads.jl
# Author: Jan van Waaij
# Created: 7 September 2023
# Last update: 7 September 2023
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

using Data, GZip, StatsBase


vcffile = "data/vcffiles_demes/two_demes_with_demes_tauC=0.25tauA=0.25numgen=1000effpop=2000_12.vcf.gz"

readfilefirstpart = "data/reads/two_demes_demes_"

numberofindividualsinthesedimentvec = 1:8 # [1, 2, 4, 8]
ϵ = 0.003
coveragevec = 8:8 # 9-20 al gedaan.

openedvcffile = GZip.open(vcffile)
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
close(openedvcffile)

for coverage in coveragevec
    for numberofindividualsinthesediment in numberofindividualsinthesedimentvec
        # I choose numberofindividualsinthesediment individuals from the ancient population. 
        individualsinthesediment = sample(ancientpopulation, numberofindividualsinthesediment, 
                                            replace=false, ordered=true)

        # I open the file where I while store the reads. 
        local openedvcffile = GZip.open(vcffile)
        # skip the header. 
        skiplines(openedvcffile, 6)
        readfile = string(readfilefirstpart, 
            numberofindividualsinthesediment, '_', coverage, "_reads_with_errors.txt.gz")
        openedreadfile = GZip.open(readfile, "w")

        comments = """
        Each line contains the number of derived reads and the total number of reads, seperated by a tab.  
        original vcf file: $vcffile
        Individuals in the sediment: $individualsinthesediment
        Total number of individuals in the sediment: $numberofindividualsinthesediment
        coverage: $coverage"""

        reads(openedvcffile, openedreadfile, individualsinthesediment, coverage, comments, ϵ)
        close(openedreadfile)
        close(openedvcffile)
    end  
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