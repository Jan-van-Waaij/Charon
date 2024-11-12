# PackageToBeShipped.jl Estimation of population parameters from environmental DNA.

PackageToBeShipped.jl is a Julia package for estimating population parameters from sediment DNA. I have tested this package extensively on all minor 1.X Julia versions (1.0-1.11). However, it is always better to use a recent version of Julia. If you have any problems, please reach out to me. 

## Install Julia (Juliaup). 
Go to [julialang.org](https://julialang.org/downloads/) and follow the installation instructions for your platform.

## Add Julia to the path. 
Check if Julia is added to your path by typing, for instance, `julia --version`. [If not, follow the instructions.](https://julialang.org/downloads/platform/)

## Install PackageToBeShipped
### Activate an environment. 
It is a good idea to start a new environment for each project [to reduce the chance of a version conflict](https://pkgdocs.julialang.org/v1/managing-packages/#conflicts). 

Press ] to activate the [package manager mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode). 
Type 

```julia
activate EnvironmentName
```

where EnvironmentName is the name you want to give to your environment. 
Next, install PackageToBeShipped to your environment: 

```julia-repl
add PackageToBeShipped
```

If you want, you can install more packages in your environment. For instance, most likely, you will need CSV and DataFrames:

```julia
] add CSV DataFrames.
```
Other valuable packages are GZip and CodecZlib for working with zipped files, StatsBase and Distributions for working with probabilities, LinearAlgebra for matrices, vectors, and related linear algebra functions, and Distributed for distributed computing.

After adding at least one package to the environment, Julia will create a new folder, `EnvironmentName`, in your current working directory. Alternatively, you could set an existing folder as the environment:
```julia
] activate path/to/the/folder
``` 
After you add at least one package, Julia will create the files Project.toml and Manifest.toml in that particular folder. Note that this action does not change the working environment. 

[Read more about Julia environments.](https://pkgdocs.julialang.org/v1/environments/)

## Data format

Our program can work with genetic data in the [DICE-2 format](https://github.com/grenaud/dice?tab=readme-ov-file#2-pop-method-input-data-format). 

A CSV file in DICE-2 format has four columns. In this order: number of ancestral reads, number of derived reads, the frequency of the anchor population, and in how many loci this combination of reads and frequency occurs. So, the first two columns contain non-negative integers, the third column is a real number between 0.0 and 1.0, and the last column is a positive integer. 

## Using the sampler. 

### Convenience script.
[Click here for an example script.](runmcmc.jl) This script sets uniform priors on `τC`, `τA`, `ϵ` and `n`, and four chains with 100'000 samples. The script saves the output as a CSV file. You can use your favourite software to analyse it. It works as follows (assuming Julia is in your path, and you use Julia 1.5 or higher):
```julia
julia --threads=4 path/to/runmcmc.jl path/to/dicefile path/to/outputfile
```
where the output is saved in `outputfile`. `dicefile` might be a CSV file or a zipped CSV file. 

### Steps to take. 
In this paragraph, I describe the steps for manually running the sampler. 
1. Save your data in the DICE-2 format as a (zipped) CSV file.
2. Start a Julia session with, for example, four threads if you want to sample four chains parallel.
```
julia --threads=4 
```
3. We must activate the environment we [created before](#activate-an-environment). You can check whether you are in the correct environment by running 
```julia
status 
```
in the package modus. If not, let's activate the environment that we created before by running
```julia
activate path/to/EnvironmentName 
```
in the package modus. Run
```julia
status
```
again to check whether the packages PackageToBeShipped, DataFrames, Distributions and CSV are installed. If any of them are lacking, please add them to your environment in the package modus via 
```julia
add PackageName
```
Alternatively, you can start your Julia session in the right environment with 
```julia 
julia --threads=4 --project=path/to/environment 
```
Make sure all the required packages are in the environment. 

4. Load the packages. 
In the REPL, run 
```julia
# load the packages. 
using PackageToBeShipped, Distributions, CSV, DataFrames
```
5. Let us specify the path to the file where our data is stored. 
```julia
dicefile = "path/to/dicefile.csv" # UNIX-systems
dicefile = "path\\to\\dicefile.csv" # Windows-sytems. 
dicefile = joinpath("path", "to", "dicefile.csv") # works on all systems. 
```
where "path/to/dicefile.csv" is the relative path to dicefile.csv from your working directory. `dicefile` might also be a zipped csv file. You can find the working directory as follows: Press ";" to come in the [shell mode](https://docs.julialang.org/en/v1/stdlib/REPL/#man-shell-mode).  
```julia-repl
shell> pwd
```
You can use `cd` to change your working directory:
```
shell> cd path/to/dicefilefolder
```
Press backspace to come back to the normal REPL. 

6. Let us specify a prior on n, (τC, τA) and ϵ. For example, you could use a uniform prior for all three:
```julia
using Distributions # Julia package for probability distributions. 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform product prior [0,1]x[0,1] on (τC, τA).
prioronn = DiscreteUniform(1, 10) # discrete uniform prior on {1, 2, ..., 10}.
prioronϵ = Uniform(0, 0.5) # uniform prior on interval[0, 0.5].
```
7. Specify how many MCMC samples you want and how many chains you want to sample. For example, if you want four chains with each 100'000 samples, you can specify 
```julia
nchains = 4 
nsteps = 100_000
```

8. We can now run our sampler, which might take a while. 
```julia
using PackageToBeShipped 
chains = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, dicefile)
```
The output of this function 
`chains` is a vector of tuples. Each tuple represents an MCMC chain. Each tuple consists of six arrays. For instance, if we want to analyse the first chain, we get 
```
nsample, τCsample, τAsample, ϵsample, accepted, logjointprob = chains[1]
```
`nsample` is the MCMC chain for the number of individuals. `τCsample`, `τAsample`, and `ϵsample` are (in this example) 100'000x10 matrices, where the k-th column is the MCMC chain of `τC`, `τA`, and `ϵ` respectively with the number of individuals conditioned on k. `accepted` and `logjointprob` are 100'000x11 matrices, where the first column of each matrix indicates whether the mcmc chain for n is accepted and the logjoint probability, respectively. The k+1-th column is the acceptance vector / log joint probability vector of the MCMC conditioned on k individuals.

You can obtain the unconditioned sample from the posterior as follows 
```julia
using PackageToBeShipped 
results = unpackposterior(chains, nsteps)
```
`results` is a DataFrame with a sample from the posterior. 

9. You can save the data frame as a CSV file as follows:
```julia
using CSV 
resultscsvfile = "path/to/results.csv" # place on your hard disk
# where you want to store your csv file. 
CSV.write(resultscsvfile, results)
```
You can use the DataFrame and the CSV file for further analysis.

## Reference 
PackageToBeShipped has three exported functions, [MCMCsampler](#mcmcsampler), [exactposterior](#exactposterior) and [unpackposterior](#unpackposterior),  that can be used after typing 
```julia
using PackageToBeShipped 
```
### MCMCsampler 
This is the MCMC sampler. It has five methods. 
```julia 
    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::AbstractVector{<:Integer}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, counts::AbstractVector{<:Integer}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::AbstractVector{<:Integer}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, df::AbstractDataFrame, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::IO, messages::Integer=nsteps÷100, scalingmessages::Bool=true, header::Bool=true)

    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::AbstractString, messages::Integer=nsteps÷100, scalingmessages::Bool=true)
```
Parameters:
* `nchains`, number of chains, a positive integer. If you want to run all your chains in parallel, start julia with number of threads equal or higher to `nchains`. 
* `nsteps` number of samples per chain, a positive integer. 
* `prioronn` the prior on n, specified as a subtype of `DiscreteUnivariateDistribution` of the Distributions Julia package. Our implementation requires that `prioronn` has support on {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, because otherwise rounding errors will accumulate too much. 
* `prioronτCτA` is the prior on (τC, τA), which allows for correlation between τC and τA. Its type is a subtype of ContinuousMultivariateDistribution. It should have support contained in [0,∞)x[0,∞). 
* `prioronϵ` is the prior on ϵ. It should have support in [0, 0.5). It is a subtype of `ContinuousUnivariateDistribution` in the Distributions package. 
* `R` a vector with coverages = ancestral reads + derived reads. Is a subtype of `AbstractVector{<:Integer}`. All coverages should be non-negative integers, and at least one should be positive. 
* `d` a vector of derived reads. Is a subtype of `AbstractVector{<:Integer}`. All elements of the vector should be non-negative integers.
* `y` a vector of frequencies. Is a subtype of AbstractVector{<:Real}. Each frequency is between 0.0 and 1.0. At least one frequency should be strictly between 0.0 and 1.0 with corresponding positive coverage. 
* `counts`, all elements should be non-negative. For each index, `counts[index]` indicates how many loci there are with `d[index]` derived reads, coverage `R[index]` and frequency `y[index]`. If `counts` is not provided, it will be calculated from `d`, `R` and `y`. 
* `messages` is an integer. If `messages` is non-positive, no message is printed. If `message` is positive, every `message` steps a message will be printed with the progress of the sampler. The default value is `nsteps÷100`. 
* `scalingmessage` is either `true` (default) or `false`. If true, a message will be printed when the scaling constant changes. 
* `df` a DataFrame from the DataFrames package in the DICE-2 format. So the first column should be the number of ancestral reads, the second column the number of derived reads, the third column the frequencies in the anchor population, and the fourth column the counts of how many times this particular combination of ancestral reads, derived reads and frequency occurs. 
* `dicefile` is either an opened DICE file, or a path to a DICE file. 
* `header` either `true` (default) or `false`. Has the dicefile a header?


The output is a vector with nchains items. Each item represents a chain. Each item is a tuple consisting of six columns, as described in step 8 of [Steps to take](#steps-to-take).

### exactposterior 

`exactposterior` is a function to calculate the posterior up to a fixed constant, only depending on the data, but not on the parameters. You can use this function for maximum likelihood estimation. `MCMCsampler` uses this function to find a good starting point for the sampler. It has two methods: 
```julia 
    exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer)
    exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, R::AbstractVector{<:Integer}, reads::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=length(n)*length(τCrange)*length(τArange)*length(ϵrange)÷100)
```
The posterior is calculated for each combination of parameters (n, τC, τA, ϵ) with n in `nrange`, τC in `τCrange`, τA in `τArange` and ϵ in `ϵrange`. So make sure that `length(nrange)*length(τCrange)*length(τArange)*length(ϵrange)` is not too large, as otherwise it will take a very long time. 

Parameters:
* `nrange` vector of n values. Subtype of AbstractVector{<:Integer}. Should be a subset of {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}. 
* `τCrange` vector of τC values. Subtype of AbstractVector{<:Real}. All values should be non-negative. 
* `τArange` vector of τA values. Subtype of AbstractVector{<:Real}. All values should be non-negative. 
* `ϵrange` vector of ϵ values. Is a subtype of AbstractVector{<:Real}. 
* `reads`

### unpackposterior 