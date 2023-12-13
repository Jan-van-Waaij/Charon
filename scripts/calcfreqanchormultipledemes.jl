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

inputfilename = "data/vcffiles_demes/multiple_demes_multipledemes_12.vcf.gz"
outputfilenames = ["data/frequencies/freqanchor/multiple_demes_freq_anchor1.txt.gz",
                   "data/frequencies/freqanchor/multiple_demes_freq_anchor2.txt.gz",
                   "data/frequencies/freqanchor/multiple_demes_freq_anchor3.txt.gz",
                   "data/frequencies/freqanchor/multiple_demes_freq_anchor4.txt.gz"]

infoblock = 9
ancientpop = infoblock .+ (1:5000)
anchorpops = [infoblock .+ (1+5000:2*5000), infoblock .+ (1+2*5000:3*5000), infoblock .+ (1+3*5000:4*5000), infoblock .+ (1+4*5000:5*5000)]

for i in 1:4
    anchorpopulation = anchorpops[i]
    outputfilename = outputfilenames[i]
    readfile = GZip.open(inputfilename)
    writefile = GZip.open(outputfilename, "w")
    # skip the header. 
    skiplines(readfile, 6);
    # save the column names. 
    comments = """
    Each line contains the frequency of the allele for the anchor population $i.
    original vcf file: $inputfilename
    output file: $outputfilename"""

    freq(readfile, writefile, anchorpopulation, comments)

    close(readfile)
    close(writefile)

    println("finished with ", writefile)
end 