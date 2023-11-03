# makediceformat.jl
# Author: Jan van Waaij
# Created: 19 December 2022
# Last updated: 16 February 2023
# Language: ≥ Julia 1.6
# 
# Purpose:
# With this script one can produce the necessary files for DICE from the reads file and the frequency file.
#
# Installation:
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open file as julia path/to/makediceformat.jl
#
# NB: works only with .vcf.gz files in the folder /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/vcffiles.
#
# Output:
# output will be saved in the folder /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/dicefiles
# the output file consist of a header with 
# Anc Der PanelFreq Num
# then rows where Num is the number of sites where there are Anc ancestrial reads, Der derived reads, PanelFreq frequency.
# The frequency is rounded to two decimals, as recommended in the DICE documentation.

using StatsBase, GZip, Data # loading packages. 

vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles") # path to the vcf files.
allfilesunfiltered = readdir(vcffilespath) # all files in the path above. 
vcffiles = String[]
for file in allfilesunfiltered # keep only the nice VCF files. 
    if file[end-6:end] != ".vcf.gz" # remove non VCF files.
        continue
    end 
    sp = split(file, '_')
    if sp[3] == "1000" # remove two_demes_1000_10_0.vcf.gz
        continue 
    end 
    push!(vcffiles, file)
end 
vcffileswithoutextension = [file[begin:end-7] for file in vcffiles] # all vcf files without extension.
readfiles = [file*"_1_20_reads.txt.gz" for file in vcffileswithoutextension] # all read files with one individual in the sediment.

readfilespath = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads" # the folder with the read files.
readfilespaths = joinpath.(readfilespath, readfiles) # all full paths to the read files. 

freqfiles = [file*"_freqanchor.txt.gz" for file in vcffileswithoutextension] # all anchor frequency files. 
freqfilespath = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqanchor/" # the path to the anchor files.  
freqfilespaths = joinpath.(freqfilespath, freqfiles) # the full paths to the frequency files. 

dicefiles = [file*"_dicefile.txt" for file in vcffileswithoutextension] # all the dice files. 
dicefilespath = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/dicefiles/" # the path to the dice files. 
dicefilespaths = joinpath.(dicefilespath, dicefiles) # all full paths to the dice files. 

for (readfile, freqfile, dicefile) in zip(readfilespaths, freqfilespaths, dicefilespaths) # make a dice file for each readfile/frequency file pair. 
    readfileio = gzopen(readfile) # open read file.
    freqfileio = gzopen(freqfile) # open frequency file. 
    derivedreadsvec, totalreadsvec = parsereadfile(readfileio) # get the derived reads and coverage from the file and save as two vectors. 
    ancientreadsvec = totalreadsvec - derivedreadsvec # get the number of ancient reads. 
    freqvec = parsefreqfile(freqfileio) # get the frequencies from the frequency file and store in a vector.
    freqvec = round.(freqvec, digits=2) # round to two decimals. 
    close(readfileio) # close the read file.
    close(freqfileio) # close the frequency file. 
    everything = Real[ancientreadsvec derivedreadsvec freqvec] # make a 3-column matrix from the three vectors. 
    # Real[] creates a Matrix{Real} type, where the first two entries are Int64's and the last entry is a Float64. 
    allrows = [everything[row, :] for row in axes(everything, 1)] # make a vector of row vectors where each element is a row from the matrix above. 
    dic = countmap(allrows) # counts the frequency of each element in allrows, saves it as a dictionary, with key the row and value the number of occurences
    # of that particular row in the matrix. 
    dicefileio = open(dicefile, "w") # open the dice file.
    print(dicefileio, "Anc\tDer\tPanelFreq\tNum") # write the header.
    for (key, value) in dic # create the dice file. key is a row of the matrix and value is the number of occurences of that particular row. 
        println(dicefileio) # new line. 
        print(dicefileio, key[1], '\t', key[2], '\t', key[3], '\t', value) # write the Anc, Der, PanelFreq, Num values. 
    end
    close(dicefileio) # close the dice file. 
end