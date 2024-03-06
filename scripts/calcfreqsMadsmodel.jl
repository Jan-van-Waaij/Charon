# calcfreqsMadsmodel.jl
# Author: Jan van Waaij
# Created: 26 February 2024
# Last update: 26 February 2024
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
# The output will have a few lines of comments, each of which preceded by #, and then each line has a number between 0.0 and 1.0, which is a frequency. 

using GZip, Data, Base.Threads

vcffile = "data/vcffilesMads/isoaltes_4pop_1kindi.vcf.gz"

populations = Dict(:CEU => 9 .+ (1:2000), :YRI => 9 .+ (2001:4000), :CHB => 9 .+ (4001:6000), :KAR => 9 .+ (6001:8000), :CEUS => 9 .+ (8001:10_000), :YRIS => 9 .+ (10_001:12_000), :CHBS => 9 .+ (12_001:14_000), :KARS => 9 .+ (14_001:16_000))

@threads for population in [:YRI, :CEU, :CHB, :KAR]
    openedvcffile = GZip.open(vcffile)

    outputfile = "data/frequencies/freqanchor/Madsmodel_freq_anchor_$(population).txt.gz"
    freqfile = GZip.open(outputfile, "w")

    anchorpopulation = populations[population]

    # skip the header. 
    skiplines(openedvcffile, 6);
    comments = """
    Each line contains the frequency of the allele for the anchor population.
    original vcf file: $vcffile
    output file: $outputfile"""

    freq(openedvcffile, freqfile, anchorpopulation, comments)

    close(openedvcffile)
    close(freqfile)
end 