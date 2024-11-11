# PackageToBeShipped.jl Estimation of population parameters from environmental DNA.

PackageToBeShipped.jl is a Julia package for estimating population parameters from sediment DNA. I have tested this package extensively, on all minor 1.X Julia versions (1.0-1.11). However, it is always better to use a recent version of Julia. For any problems, please reach out to me. 

## Install Julia (Juliaup). 
Go to [julialang.org](https://julialang.org/downloads/) and follow the installation instructions for your platform.

## Add Julia to path. 
Check if Julia is added to your path by typing, for instance `julia --version`.  [If not, follow the instructions.](https://julialang.org/downloads/platform/)

## Install PackageToBeShipped
### Activate an environment. 
For beginners of Julia, it might be a good idea to 

Press ] to activate the [package manager mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode). 
Type 

```julia
activate EnvironmentName
```

where EnvironmentName is the name you want to give to your environment. 
Next install PackageToBeShipped to your environment: 

```julia
add PackageToBeShipped
```

If you want, you can install more packages to your environment. For instance, most likely, you will need CSV, DataFrames,
GZip, CodecZlib, StatsBase, LinearAlgebra, Distributions and Distributed.

```julia
] add CSV DataFrames GZip CodecZlib StatsBase LinearAlgebra Distributions Distributed
```

After you added at least one package to the environment, Julia will create a new folder `EnvironmentName` in your current working directory. Alternatively, you could set an existing folder as the environment:
```julia
] activate path/to/the/folder
``` 
After you add at least one package, Julia will create the files Project.toml and Manifest.toml in that particular folder. Note that this action does not change the working environment. 

[Read more about Julia environments.](https://pkgdocs.julialang.org/v1/environments/)

## Data format. 

Our program can work with genetic data in the [DICE-2 format](https://github.com/grenaud/dice?tab=readme-ov-file#2-pop-method-input-data-format). 

This is a CSV-file with four columns. In this order: number of ancestral reads, number of derived reads, the frequency of the anchor population, and in how many loci this combination of reads and frequency occurs. So the first two columns contain non-negative integers, the third column a real number between 0.0 and 1.0 and the last column is a positive integer. 

## Using the sampler. 

### Convenience script.
[Click here for an example script.](runmcmc.jl) This script sets uniform priors on `τC`, `τA`, `ϵ` and `n`, and 4 chains with 100'000 samples. The output is saved as a CSV file. You can use your favourite software to analyse it. It works as follows (assuming Julia is in your path, and you use Julia 1.5 or higher):
```julia
julia --threads=4 path/to/runmcmc.jl path/to/dicefile path/to/outputfile
```
where the output is saved in `outputfile`. `dicefile` might be a CSV file or a zipped CSV file. 

### Steps to take. 
This are the steps that you need to take. 
1. Put your data in the DICE-2 format, and save it as a (zipped) CSV file. 
2. Specify a prior on n, (τC, τA) and ϵ. For example you could use a uniform prior for all three:
```julia
using Distributions # Julia package for probability distributions. 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform product prior [0,1]x[0,1] on (τC, τA).
prioronn = DiscreteUniform(1, 10) # discrete uniform prior on {1, 2, ..., 10}.
prioronϵ = Uniform(0, 0.5) # uniform prior on interval[0, 0.5].
```
3. Specify how many mcmc samples you want and how many chains you want to sample. For example, if you want 4 chains with each 100'000 samples you can specify 
```julia
nchains = 4 
nsteps = 100_000
```
4. Next you need to set the path to your DICE file. For example 
```julia
dicefile = "path/to/dicefile.csv" # UNIX-systems
dicefile = "path\\to\\dicefile.csv" # Windows-sytems. 
```
where "path/to/dicefile.csv" is the relative path to dicefile.csv from your working directory. You can find the working directory as follows: Press ";" to come in the [shell mode](https://docs.julialang.org/en/v1/stdlib/REPL/#man-shell-mode).  
```julia-repl
shell> pwd
shell> cd path/to/dicefilefolder
```
Press backspace to come back in the normal REPL. 

5. If you want to save the MCMC results as a CSV file, set a path to where you want to store the results. Make sure it is not an existing file, in case you don't want the data to be overwritten. 
```julia
results = "path/to/results.csv"
```

6. Congratulations! Now we can run our code. 
Let's activate the environment that we created before:
Press ] and type the following in the package REPL:
```julia
activate path/to/EnvironmentName 
```
Type 
```julia
status
```
to check whether the packages PackageToBeShipped, Distributions and CSV are installed. If one lacks, please add it in the package-REPL via 
```julia
pkg> add PackageName
```
in the Package-REPL. 
Press escape to come in the REPL and copy-paste the following code: 
```julia
using PackageToBeShipped, Distributions, CSV, DataFrames

prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]

nchains = 4 
nsteps = 100_000

dicefile = "path/to/dicefile.csv"
resultscsvfile = "path/to/results.csv"

chains = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, dicefile)
```

`chains` is a vector of tuples. Each tuple represents a MCMC-chain. Each tuple consists of six arrays. For instance if we want to analyse the first chain, we get 
```
nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = chains[1]
```
`nsample` is the MCMC chain for the number of individuals. `τCsample`, `τAsample`, and `ϵsample` are (in this example) 100'000x10 matrices, where the k-th column is MCMC chain of `τC`, `τA`, and `ϵ` respectively with number of individuals conditioned on k. `acceptedvec` and `logjointprob` are 100'000x11 matrices, where the first column of each matrix indicates whether the mcmc chain for n is accepted and the logjoint probability. The k+1-th column is the acceptance vector / log joint probability vector of the mcmc conditioned on the number of k.

You can distract a sample of the posterior as follows: 
```julia
results = DataFrame(nsample = Int[],
            τCsample = Float64[],
            τAsample = Float64[],
            ϵsample = Float64[],
            accepted = Bool[],
            logjointprob = Float64[],
            chainid = Int[]
            )
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
    df = DataFrame(nsample = nsample,
                        τCsample = τCsampleunpacked,
                        τAsample = τAsampleunpacked,
                        ϵsample = ϵsampleunpacked,
                        accepted = acceptedvecunpacked,
                        logjointprob = logjointprob[:,1],
                        chainid = fill(Int(index), nsteps)
                        )
    results = vcat(results, df)
end 
```



Congratulations! You got your a sample from your posterior. You can save the dataframe as CSV file as follows:
```julia
CSV.write(resultscsvfile, results)
```

 They are stored in the Julia session as a DataFrame and on the disk as a CSV file, which you can use for further analysis. 
