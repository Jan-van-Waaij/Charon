using Pkg 


Pkg.activate("Runmcmc")
#Pkg.add(["PackageToBeShipped", "Distributions", "CSV", "DataFrames"])
Pkg.add(["Distributions", "CSV", "DataFrames"])

using PackageToBeShipped, Distributions, CSV, DataFrames 

if length(ARGS) != 2 
    println("You should provide two arguments: (1) path to the dice file and (2) path to the results.")
    exit(1)
end 

dicefile = ARGS[1]
resultsfile = ARGS[2]

if !isfile(dicefile)
    println("No file found at $dicefile.")
    exit(2)
end 

# specify the number of chains and the number of samples. 
nchains = 4 
nsteps = 100_000

# specify the prior. 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform product prior [0,1]x[0,1] on (τC, τA).
prioronn = DiscreteUniform(1, 10) # discrete uniform prior on {1, 2, ..., 10}.
prioronϵ = Uniform(0, 0.5) # uniform prior on interval[0, 0.5].

chains = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, dicefile)

println(typeof(chains))

df = DataFrame(nsample = Int[],
            τCsample = Float64[],
            τAsample = Float64[],
            ϵsample = Float64[],
            accepted = Bool[],
            logjointprob = Float64[],
            chainid = Int[]
            )
CSV.write(resultsfile, df)

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
    local df = DataFrame(nsample = nsample,
                        τCsample = τCsampleunpacked,
                        τAsample = τAsampleunpacked,
                        ϵsample = ϵsampleunpacked,
                        accepted = acceptedvecunpacked,
                        logjointprob = logjointprob[:,1],
                        chainid = fill(Int(index), nsteps)
                        )
    CSV.write(resultsfile, df, append=true)
end 
