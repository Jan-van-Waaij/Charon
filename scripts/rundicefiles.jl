# This is a file that users can use to execute a dicefile. 

# Steps:
# install Julia
# add Julia to path. 
# Run this script as follows
# julia path/to/rundicefiles.jl path/to/dicefile.txt nchains nsamples path/to/outputfile.txt

# dicefile is the dicefile.  
# nchains is the number of chains
# nsamples is the number of mcmc steps. 
# outputfile will contain the mcmc output. 

# DICE 2-pop format:
# ancestral_reads derived_reads frequency counts 
# See https://github.com/grenaud/dice?tab=readme-ov-file#2-pop-method-input-data-format


using Pkg
Pkg.activate(".")
Pkg.add(["StatsBase", "DataFrames", "CSV", "GZip", "PosteriorSampling"])

using StatsBase, DataFrames, CSV, GZip, PosteriorSampling

# check the arguments. 
if length(ARGS) != 4
    println("You should give 4 arguments: dicefile, nchains, nsteps, outputfile")
    exit(1)
end 
dicefile = ARGS[1]
if !isfile(dicefile)
    println("$dicefile does not exist.")
    exit(2)
end 
if !all(isdigit, ARGS[2])
    println("The second argument nchains should be an integer")
    exit(3)
end 
if !all(isdigit, ARGS[3])
    println("The second argument nsteps should be an integer")
    exit(4)
end 

nchains = parse(Int, ARGS[2])
if nchains ≤ 0 
    println("The second argument nchains should be a positive number")
    exit(5)
end
nsamples = parse(Int, ARGS[3])
if nsample ≤ 0 
    println("The third argument nsteps should be a positive number")
    exit(6)
end 
outputfile = ARGS[4]

prioronn = DiscreteUniform(1,10)
prioronτCτA = Product([Uniform(0,1), Uniform(0,1)])
priononϵ = Uniform(0,0.5)
results = MCMCSimProp(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, dicefile)

CSV.write(outputfile, results)

