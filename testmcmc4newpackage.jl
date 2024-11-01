using PackageToBeShipped 
using Data 
using DataFrames
using CSV 
using Distributions 
using Dates 
using StatsBase 
using GZip 


ancientpop = "CHBS"
anchorpop = "YRI"
truenumind = 2 
coverage = 5 

trueϵ = 0.003

nsteps = 100_000 # number of steps in mcmc chain 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]
nchains = 4 # number of chains. 

# test 1 
#= 
csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "MadsModel", "test.csv")
zippedcsvfile = csvfile * ".gz"

starttime = now()
frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", "freqanchor", "Madsmodel_freq_anchor_$(anchorpop).txt.gz")
freqs = parsefreqfile(frequencyfile)
readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", string("Mads_model_pop_", ancientpop, "_num_ind_",
truenumind, "_cov_", coverage, "_reads_with_errors.txt.gz"))
derivedreadsvec, coveragevec = parsereadfile(readsfile)
@assert length(unique(coveragevec)) == 1
@assert coveragevec[1] == coverage
@assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

df = DataFrame(anchorpop = String[],
            ancientpop = String[],
            trueϵ = Float64[],
            truenumind = Int[],
            coverage = Int[],
            nsample = Int[],
            τCsample = Float64[],
            τAsample = Float64[],
            ϵsample = Float64[],
            accepted = Bool[],
            logjointprob = Float64[],
            chainid = Int[]
            )
CSV.write(csvfile, df)
#results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
chains = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coveragevec, derivedreadsvec, freqs, 0, false)
for index in eachindex(chains)
    nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = chains[index]
    τCsampleunpacked = Vector{Float64}(undef, nsteps)
    τAsampleunpacked = Vector{Float64}(undef, nsteps)
    ϵsampleunpacked = Vector{Float64}(undef, nsteps)
    for index in 1:nsteps
        τCsampleunpacked[index] = τCsample[index, nsample[index]]
        τAsampleunpacked[index] = τAsample[index, nsample[index]]
        ϵsampleunpacked[index] = ϵsample[index, nsample[index]]
    end 
    acceptedvecunpacked = BitVector(undef, nsteps)
    acceptedvecunpacked[1] = true 
    for index in 2:nsteps
        everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) && (nsample[index] == nsample[index-1])
        acceptedvecunpacked[index] = !everythingthesame
    end 
    local df = DataFrame(
                        anchorpop = fill(anchorpop, nsteps),
                        ancientpop = fill(ancientpop, nsteps),
                        trueϵ = fill(trueϵ, nsteps),
                        truenumind = fill(truenumind, nsteps),
                        coverage = fill(coverage, nsteps),
                        nsample = nsample,
                        τCsample = τCsampleunpacked,
                        τAsample = τAsampleunpacked,
                        ϵsample = ϵsampleunpacked,
                        accepted = acceptedvecunpacked,
                        logjointprob = logjointprob[:,1],
                        chainid = fill(Int(index), nsteps)
                        )
    CSV.write(csvfile, df, append=true)
end  
command = `gzip -f $csvfile`
run(command)
finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("It took ", numhours, " and ", remainingtimeminutes, '.')
println("Finished with ", csvfile)
=#

# test 2 
#= 
csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "MadsModel", "test2.csv")
zippedcsvfile = csvfile * ".gz"

starttime = now()
frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", "freqanchor", "Madsmodel_freq_anchor_$(anchorpop).txt.gz")
freqs = parsefreqfile(frequencyfile)
readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", string("Mads_model_pop_", ancientpop, "_num_ind_",
truenumind, "_cov_", coverage, "_reads_with_errors.txt.gz"))
derivedreadsvec, coveragevec = parsereadfile(readsfile)
@assert length(unique(coveragevec)) == 1
@assert coveragevec[1] == coverage
@assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

# make a Dice file. 
ancestralreads = coveragevec - derivedreadsvec 
cm = countmap(zip(ancestralreads, derivedreadsvec, freqs))
ancestralreads = [key[1] for key in keys(cm)]
derivedreadsvec = [key[2] for key in keys(cm)]
freqs = [key[3] for key in keys(cm)]
counts = collect(values(cm))
dicedf = DataFrame(ancestralreads = ancestralreads, derivedreads = derivedreadsvec, freq = freqs, count = counts)


df = DataFrame(anchorpop = String[],
            ancientpop = String[],
            trueϵ = Float64[],
            truenumind = Int[],
            coverage = Int[],
            nsample = Int[],
            τCsample = Float64[],
            τAsample = Float64[],
            ϵsample = Float64[],
            accepted = Bool[],
            logjointprob = Float64[],
            chainid = Int[]
            )
CSV.write(csvfile, df)
#results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
chains = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, dicedf, 0, false)
for index in eachindex(chains)
    nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = chains[index]
    τCsampleunpacked = Vector{Float64}(undef, nsteps)
    τAsampleunpacked = Vector{Float64}(undef, nsteps)
    ϵsampleunpacked = Vector{Float64}(undef, nsteps)
    for index in 1:nsteps
        τCsampleunpacked[index] = τCsample[index, nsample[index]]
        τAsampleunpacked[index] = τAsample[index, nsample[index]]
        ϵsampleunpacked[index] = ϵsample[index, nsample[index]]
    end 
    acceptedvecunpacked = BitVector(undef, nsteps)
    acceptedvecunpacked[1] = true 
    for index in 2:nsteps
        everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) && (nsample[index] == nsample[index-1])
        acceptedvecunpacked[index] = !everythingthesame
    end 
    local df = DataFrame(
                        anchorpop = fill(anchorpop, nsteps),
                        ancientpop = fill(ancientpop, nsteps),
                        trueϵ = fill(trueϵ, nsteps),
                        truenumind = fill(truenumind, nsteps),
                        coverage = fill(coverage, nsteps),
                        nsample = nsample,
                        τCsample = τCsampleunpacked,
                        τAsample = τAsampleunpacked,
                        ϵsample = ϵsampleunpacked,
                        accepted = acceptedvecunpacked,
                        logjointprob = logjointprob[:,1],
                        chainid = fill(Int(index), nsteps)
                        )
    CSV.write(csvfile, df, append=true)
end  
command = `gzip -f $csvfile`
run(command)
finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("It took ", numhours, " and ", remainingtimeminutes, '.')
println("Finished with ", csvfile)
=#

# test 3 

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "MadsModel", "test2.csv")
zippedcsvfile = csvfile * ".gz"

starttime = now()
frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", "freqanchor", "Madsmodel_freq_anchor_$(anchorpop).txt.gz")
freqs = parsefreqfile(frequencyfile)
readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", string("Mads_model_pop_", ancientpop, "_num_ind_",
truenumind, "_cov_", coverage, "_reads_with_errors.txt.gz"))
derivedreadsvec, coveragevec = parsereadfile(readsfile)
@assert length(unique(coveragevec)) == 1
@assert coveragevec[1] == coverage
@assert length(freqs) == length(derivedreadsvec) == length(coveragevec)

# make a Dice file. 
ancestralreads = coveragevec - derivedreadsvec 
cm = countmap(zip(ancestralreads, derivedreadsvec, freqs))
ancestralreads = [key[1] for key in keys(cm)]
derivedreadsvec = [key[2] for key in keys(cm)]
freqs = [key[3] for key in keys(cm)]
counts = collect(values(cm))
dicedf = DataFrame(ancestralreads = ancestralreads, derivedreads = derivedreadsvec, freq = freqs, count = counts)
CSV.write("testdicefile.txt", dicedf)

file = open("testdicefile.txt")


df = DataFrame(anchorpop = String[],
            ancientpop = String[],
            trueϵ = Float64[],
            truenumind = Int[],
            coverage = Int[],
            nsample = Int[],
            τCsample = Float64[],
            τAsample = Float64[],
            ϵsample = Float64[],
            accepted = Bool[],
            logjointprob = Float64[],
            chainid = Int[]
            )
CSV.write(csvfile, df)
#results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, coverage, derivedreadsvec, freqs, 2000)
chains = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, file, 0, false)
for index in eachindex(chains)
    nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = chains[index]
    τCsampleunpacked = Vector{Float64}(undef, nsteps)
    τAsampleunpacked = Vector{Float64}(undef, nsteps)
    ϵsampleunpacked = Vector{Float64}(undef, nsteps)
    for index in 1:nsteps
        τCsampleunpacked[index] = τCsample[index, nsample[index]]
        τAsampleunpacked[index] = τAsample[index, nsample[index]]
        ϵsampleunpacked[index] = ϵsample[index, nsample[index]]
    end 
    acceptedvecunpacked = BitVector(undef, nsteps)
    acceptedvecunpacked[1] = true 
    for index in 2:nsteps
        everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) && (nsample[index] == nsample[index-1])
        acceptedvecunpacked[index] = !everythingthesame
    end 
    local df = DataFrame(
                        anchorpop = fill(anchorpop, nsteps),
                        ancientpop = fill(ancientpop, nsteps),
                        trueϵ = fill(trueϵ, nsteps),
                        truenumind = fill(truenumind, nsteps),
                        coverage = fill(coverage, nsteps),
                        nsample = nsample,
                        τCsample = τCsampleunpacked,
                        τAsample = τAsampleunpacked,
                        ϵsample = ϵsampleunpacked,
                        accepted = acceptedvecunpacked,
                        logjointprob = logjointprob[:,1],
                        chainid = fill(Int(index), nsteps)
                        )
    CSV.write(csvfile, df, append=true)
end  
command = `gzip -f $csvfile`
run(command)
finishtime = now()
timedifference = finishtime - starttime
numhours = round(timedifference, Hour, RoundDown)
remainingtime = timedifference - numhours
remainingtimeinmilliseconds = Dates.Millisecond(Dates.toms(remainingtime))
remainingtimeminutes = round(remainingtimeinmilliseconds, Minute)
println("It took ", numhours, " and ", remainingtimeminutes, '.')
println("Finished with ", csvfile)

close(file)