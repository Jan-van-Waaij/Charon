"""
    module Data 

This package contains functions for data simulation and processing: creating simulated read data from VCF files, calculating frequencies from VCF files, and store them in files, along with functions to parse read and frequency files.
"""
module Data

# Write test functions. 

# Required packages. 
using Distributions, # for distributions. Here we use the Binomial distribution and the rand function.
      DelimitedFiles, # To work with tab seperated files. Here we use the functions readdlm and writedlm.
      GZip # To open zipped vcf files. 

import Base: showerror # format the error and print it. 

# Exported functions. Available once you add "using Data" to your code. 
# Other functions can be used with Data.functionname(arguments), which requires at least "import Data" or "using Data". 
export skiplines, # skip the first n lines of a file
       reads, # generate simulated read data from a vcf file. 
       freq, # calculate the frequency of a certain allele of a population from a vcf file.
       parsefreqfile, # parse the frequency file
       parsereadfile, # parse the read file 
       addhashtagstoeachline # add '#' to each line in the string. 

# Define custom DataError Error exception.
# test written.
"""
    DataError

DataError exception. Something is wrong with the data format. 
"""
struct DataError <: Exception
    message::String
end
# Constructor for DataError.
"""
    DataError(x...)

Data is in the wrong format, details are specified in x. 
"""
DataError(x...) = DataError(string(x...)) 
# Formatting of the error message. 
Base.showerror(io::IO, e::DataError) = print(io, "DataError: ", e.message)

# Define custom FileError Error exception.
# test written.
struct FileError <: Exception
    message::String
end
# Constructor for FileError.
"""
    FileError(x...)

Somethings is wrong with the file.
"""
FileError(x...) = FileError(string(x...)) 
# Formatting of the error message. 
Base.showerror(io::IO, e::FileError) = print(io, "FileError: ", e.message)

# Try to find a library function that does this.
# test written.
"""
    skiplines(openedfile::IO, numberoflinestobeskipped::Int)

Skip the first numberoflinestobeskipped lines of openedfile. 
"""
function skiplines(openedfile::IO, numberoflinestobeskipped::Int)
    numberoflinestobeskipped >= 0 || throw(DomainError("numberoflinestobeskipped must be nonnegative"))
    for _ in Base.OneTo(numberoflinestobeskipped)
        readline(openedfile)
    end
end

# test written.
"""
    isderived(x::Char)

Assumption: x is a character in "0|0", "0|1", "1|0", or "1|1". Returns true when x == '1', and zero otherwise. Throws a DataError when x is not '0', '1' or '|'.  
"""
function isderived(x::Char)
    if x == '1'
        return true
    elseif x == '0'
        return false 
    elseif x == '|'
        return false 
    else 
        throw(DataError("Expected '0', '1', or '|', but got ", x, '.'))
    end 
end 

# test written
"""
    countderivedalleles(line::Vector{SubString{String}}, population::AbstractVector{Int})

Our data comes in the form of a vector of substrings, where each substring is either "0|0", "0|1", "1|0" or "1|1", 
which adds 0, 1, 1 or 2 derived alleles to the total number of derived alleles, respectively. This function counts the total number of derived alleles of the population at the specific locus to which the line corresponds.
"""
countderivedalleles(line::Vector{SubString{String}}, population::AbstractVector{Int}) = sum(count.(isderived, line[population]))
# line[population] selects all genotypes belonging to individuals in the population. 
# count is broadcasted (indicated by the . after the function name) over all SubStrings in line[population]. 
# Each substring is either "0|0", "0|1", "1|0", or "1|1". The function 'count' counts how many times each substring
# has characted '1' (so 0, 1, 1, or 2 times respectively). Next, the results for each individual are summed,
# over all individuals in the population, which gives the total number of derived alleles in the population.

# test written
"""
    addhashtagstoeachline(str::String)

Precedes each line in str with #, except for an empty string, which is left untouched. 
"""
function addhashtagstoeachline(str::String)
    if str == ""
        return "" # leave an empty string untouched. 
    else
        sp = split(str, '\n') # split the string by newline.
        vec = '#'.*sp.*'\n' # add # before and \n after each line. 
        return string(vec...) # concatenate the vector of strings to one string and return. 
    end     
end 

# tested.
"""
    openfile(file::String)
    openfile(file::String, fileextension::String)

Open file in read modus, checks whether it exists and, optionally, if it ends with fileextension. 
"""
function openfile(file::String)
    #endswith(file, fileextension) || throw(FileError("$file is not a $fileextension file."))
    isfile(file) || throw(FileError("$file does not exist."))
    return GZip.open(file) 
end
function openfile(file::String, fileextension::String)
    endswith(file, fileextension) || throw(FileError("$file is not a $fileextension file."))
    return openfile(file)
end

"""
    addcomments(file::IO, comments::String)

Add comments to the file, and precede them by #'s. Do nothing when comments is the empty string.
"""
function addcomments(file::IO, comments::String)
    if comments != ""
        commentswithhashtags = addhashtagstoeachline(comments)
        print(file, commentswithhashtags)
    end 
end

# tests written. 
"""
    openoutputfileandwritecomments(outputfile::String, overwrite::Bool, comments::String)

Open outputfile and write the comments, next return the IO stream. Throws a FileError when overwrite is false, which prevents accidentally overwriting files. It checks whether writefile is a gzipped (.gz) file.
"""
function openoutputfileandwritecomments(outputfile::String, overwrite::Bool, comments::String)
    isfile(outputfile) && !overwrite && throw(FileError("$outputfile does already exist."))
    endswith(outputfile, ".gz") || throw(FileError("$outputfile should be a zip file (so with file extension .gz)."))
    openedoutputfile = GZip.open(outputfile, "w")
    addcomments(openedoutputfile, comments)
    return openedoutputfile
end

# tests written 
"""
    reads(numberofderivedalles::Int, numberofindividuals::Int, coverage::Int, error::Float64=0.0)

Simulates the number of derived reads, with coverage and numberofderivedalles and numberofindividuals. It is binomially distributed with parameters 
coverage and numberofderivedalles / (2*numberofindividuals). Note that we assume diploid individuals. With probability error, a derived read is changed into an ancestral read or vice versa an ancestral read is changed into a derived read.
"""
function reads(numberofderivedalles::Int, numberofindividuals::Int, coverage::Int, error::Float64=0.0)
    coverage ≥ 0 || throw(ArgumentError("Coverage should be a non-nonegative integer, you provided $coverage."))
    0.0 ≤ error ≤ 1.0 || throw(ArgumentError(string("The error probability should be between 0.0 and 1.0, but you provided ", error, ".")))
    distderivedreads = Binomial(coverage, numberofderivedalles / (2*numberofindividuals))
    numderivedreads = rand(distderivedreads)
    # derived reads remains derived; no error occurs. 
    disterrorderivedreads = Binomial(numderivedreads, 1-error)
    # ancestral read becomes derived through error.
    disterrorancestralreads = Binomial(coverage-numderivedreads, error)
    numderivedreads = rand(disterrorderivedreads) + rand(disterrorancestralreads)
    return numderivedreads, coverage
end

# tests written.
"""
    reads(numberofderivedalles::Int, numberofindividuals::Int, distcoverage::Distribution{Univariate, Discrete}, error::Float64=0.0)

Simulate the number of derived reads and total number of reads. Returns a tuple with number of derived reads and total number of reads. 
distnumreads is the distribution of the number of reads, for instance Poisson(10). Returns a tuple of Ints where 
the first argument is the number of derived reads and the second number is the total number of reads. With probability error, a derived read is changed into an ancestral read or vice versa an ancestral read is changed into a derived read.
"""
function reads(numberofderivedalles::Int, numberofindividuals::Int, distcoverage::Distribution{Univariate, Discrete}, error::Float64=0.0)
    coverage = rand(distcoverage)
    return reads(numberofderivedalles, numberofindividuals, coverage, error)
end

# test written. 
"""
    reads(vcffile::IO, individualsinthesediment::Vector{Int}, coverage, error::Float64=0.0)

vcffile is an opened (zipped) VCF file (.vcf or .vcf.gz). Select the individuals at positions individualsinthesediment and count the number of derived alleles of all individuals together, at each locus. Then simulate reads. The parameter coverage is either an Int or a Distribution{Univariate, Discrete} distribution.

Output: two vectors (Vector{Int}) where the first one are the number of derived reads for each locus, and the second is the coverage for each locus. In case coverage is an Int, then the second vector is constant. 

With probability error, a derived read is changed into an ancestral read or vice versa an ancestral read is changed into a derived read.
"""
function reads(vcffile::IO, individualsinthesediment::Vector{Int}, coverage, error::Float64=0.0)
    allderivedreads = Int[] # vector with number of derived reads.
    allcoverages = Int[] # vector with coverages. 
    numberofindividualsinthesediment = length(individualsinthesediment)
    line = readline(vcffile) # read the first line of the vcf file.
    # each line is a line of genotypes at a certain locus. 
    while line != "" # continue until you have read the whole file. 
        sp = split(line) # split the line by tabs. 
        derivedallelecountssediment = countderivedalleles(sp, individualsinthesediment) 
        # the substrings at the positions in individualsinthesediment are of the form "0|0", "0|1", "1|0", or "1|1".
        numberofderivedreads, simcoverage = reads(derivedallelecountssediment, numberofindividualsinthesediment, coverage, error)
        # simulate reads and coverage. 
        push!(allderivedreads, numberofderivedreads) # push the results to the vectors. 
        push!(allcoverages, simcoverage)
        line = readline(vcffile) 
    end 
    return allderivedreads, allcoverages
end  

# test written 
"""
    reads(vcffile::String, individualsinthesediment::Vector{Int}, coverage, error::Float64=0.0)

Open VCF file vcffile and generate reads. Returns two vectors, the first containing the total number of derived reads and the second containing all coverages. This function closes the VCF file afterwards. 

With probability error, a derived read is changed into an ancestral read or vice versa an ancestral read is changed into a derived read.
"""
function reads(vcffile::String, individualsinthesediment::Vector{Int}, coverage, error::Float64=0.0)
    openedvcffile = openfile(vcffile, ".vcf.gz")
    allderivedreads, allcoverages = reads(openedvcffile, individualsinthesediment, coverage, error)
    close(openedvcffile)
    return allderivedreads, allcoverages
end 

# test written 
"""
    reads(vcffile::IO, readsfile::IO, individualsinthesediment::Vector{Int}, coverage, comments::String="", error::Float64=0.0)

Simulate reads with reads(readfile, individualsinthesediment, coverage) and write them to writefile. writefile contains two columns, the first containing the total number of derived reads and the second containing all coverages.  The comments (if any) are preceded with '#' and written to readsfile.

With probability error, a derived read is changed into an ancestral read or vice versa an ancestral read is changed into a derived read.
"""
function reads(vcffile::IO, readsfile::IO, individualsinthesediment::Vector{Int}, coverage, comments::String="", error::Float64=0.0)
    addcomments(readsfile, comments)
    allderivedreads, allcoverages = reads(vcffile, individualsinthesediment, coverage, error) # simulate reads.
    writedlm(readsfile, [allderivedreads allcoverages]) # write the results to writefiles.
end

# test written 
"""
    reads(vcffile::String, readsfile::String, individualsinthesediment::Vector{Int}, coverage, overwrite::Bool=false, comments::String="", error::Float64=0.0)

Simulate reads with reads(readfile, individualsinthesediment, coverage) and write them to writefile. writefile contains two columns, the first containing the total number of derived reads and the second containing all coverages. vcffile should be an vcf file with .vcf.gz extension and readsfile should be a zip file with .gz extension. The comments (if any) are preceded with '#' and written to readsfile. overwrite decides whether readsfile might be overwritten in case it already exist. 

With probability error, a derived read is changed into an ancestral read or vice versa an ancestral read is changed into a derived read.
"""
function reads(vcffile::String, readsfile::String, individualsinthesediment::Vector{Int}, coverage, overwrite::Bool=false, comments::String="", error::Float64=0.0)
    openedvcffile = openfile(vcffile, ".vcf.gz")
    openedreadsfile = openwritefileandwritecomments(readsfile, overwrite, comments)
    allderivedreads, allcoverages = reads(openedvcffile, individualsinthesediment, coverage, error)
    writedlm(openedreadsfile, [allderivedreads allcoverages])
    close(openedvcffile)
    close(openedreadsfile)
end

# test written 
# This function is needed so that the compiler can distinghuish between all the options. 
#reads(vcffile::String, readsfile::String, individualsinthesediment::Vector{Int}, coverage, comments::String="", error::Float64=0.0) =
#    reads(vcffile, readsfile, individualsinthesediment, coverage, false, comments, error)

#= HOLA PROBLEEMPJE!!!!!
IK KREEG DE VOLGENDE FOUTMELDING
WARNING: Method definition reads(String, String, Array{Int64, 1}, Any) in module Data at /home/people/jvawa/.julia/packages/Data/ibwf7/src/Data.jl:253 overwritten at /home/people/jvawa/.julia/packages/Data/ibwf7/src/Data.jl:264.

OPLOSSING 1: REMOVER 264 EN KIJK WAT ER GEBEURT!
=#

#=
Algemene stramien: 
voeg # toe aan comments. 
schrijf de comments. 

initialise some variables.

readline 

loop over all lines and write all the results to some other file, or store them into a vector. 
=#

# tests written
"""
    calc_freq(line::String, population::UnitRange{Int})

Caculate the frequency of the derived allele at a certain locus, by counting the number of derived alleles and dividing it by twice the population size (because we have diploid individuals). 
"""
function calc_freq(line::String, population::UnitRange{Int})
    sp = split(line) # split the string by tabs. 
    numalleles = 2*length(population) # we have diploid individuals, so the number of alleles is twice the 
    # number of individuals. 
    numderalleles = countderivedalleles(sp, population) 
    return numderalleles / numalleles
end

# test written 
"""
    freq(openedzipfile::IO, population::UnitRange{Int})

The file is an opened zipped VCF file (.vcf.gz). For each locus the frequency of the derived allele is calculated. The function returns a vector with all frequencies, where the k-th element corresponds to the frequency of the k-th locus.  
"""
function freq(openedzipfile::IO, population::UnitRange{Int})
    allfreqs = Float64[] # vector with all frequencies.
    line = readline(openedzipfile) # read the first line of the vcf file.
    while line != "" # continue untill the file is finished. 
        fr = calc_freq(line, population) # caculate the frequency at this locus in this population.
        push!(allfreqs, fr) # add the result to the vector. 
        line = readline(openedzipfile) # read the next line. 
    end 
    return allfreqs
end

# test written 
"""
    freq(vcffile::String, population::UnitRange{Int})

The file is a zipped VCF file (.vcf.gz). For each locus the frequency of the derived allele is calculated. 
The function returns a vector with all frequencies, where the k-th element corresponds to the frequency of the k-th locus.  
"""
function freq(vcffile::String, population::UnitRange{Int})
    openedvcffile = openfile(vcffile, ".vcf.gz")
    allfreqs = freq(openedvcffile, population)
    close(openedvcffile)
    return allfreqs
end

# test written 
"""
    freq(readfile::IO, writefile::IO, population::UnitRange{Int}, comments::String="")

Apply freq(openedzipfile, population) and write the results to writefile.
"""
function freq(readfile::IO, writefile::IO, population::UnitRange{Int}, comments::String="")
    commentswithhashtags = addhashtagstoeachline(comments) # precede each line of the comments with a '#'. 
    println(writefile, commentswithhashtags) # write the comments to the writefile.
    allfreqs = freq(readfile, population) # calculate all the frequencies in the readfile. 
    writedlm(writefile, allfreqs) # write the results to the writefile.
end

# test written 
function freq(vcffile::String, freqfile::String, population::UnitRange{Int}, overwrite::Bool=false, comments::String="")
    allfreqs = freq(vcffile, population)
    openedfreqfile = openwritefileandwritecomments(freqfile, overwrite, comments)
    writedlm(openedfreqfile, allfreqs) # write the results to the freqfile.
    close(openedfreqfile)
end

# test written 
# additional method for when the second last argument 
# (with default argument false) is omited, but the last  
# method not. Gives MethodError when omited.  
freq(vcffile::String, freqfile::String, population::UnitRange{Int}, comments::String) = 
    freq(vcffile, freqfile, population, false, comments)

# tests written  
"""
    parsefreqfile(readfile::IO)

Parse a frequency file and return a vector with frequencies. 
"""
parsefreqfile(readfile::IO) = vec(readdlm(readfile, '\t', Float64, comments = true))

# tests written 
"""
    parsefreqfile(readsfile::String)

Parse the frequency file and return a vector with frequencies. readsfiles is a .txt.gz file. 
"""
function parsefreqfile(freqsfile::String)
    openedfreqsfile = openfile(freqsfile)#, ".txt.gz")
    allfreqs = parsefreqfile(openedfreqsfile)
    close(openedfreqsfile)
    return allfreqs
end

# Test written.
"""
    parsereadfile(readfile::IO)

Parse the read files and return two vectors, the first one containing the number of derived reads for each locus, and the second containing the coverage for each locus. 
"""
function parsereadfile(readfile::IO, delim = '\t')
    M = readdlm(readfile, delim, Int, comments = true)
    return M[:,1], M[:,2]
end 

# Test written.
"""
    parsereadfile(readsfile::String)

Parse the reads file and return two vectors, the first containing the number of derived reads and the second the coverages. The file is assumed to be a .txt.gz file. 
"""
function parsereadfile(readsfile::String, delim = '\t')
    openedreadsfile = openfile(readsfile)#, ".txt.gz")
    derivedreads, coverages = parsereadfile(openedreadsfile, delim)
    close(openedreadsfile)
    return derivedreads, coverages
end

end # module Data
