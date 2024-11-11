# PackageToBeShipped.jl Estimation of population parameters from environmental DNA.

## Install Julia (Juliaup). 
Go to [julialang.org](julialang.org) and follow the installation instructions there for your platform. 

## Add Julia to path. 
Check if Julia is added to your path by typing, for instance `julia --version`. If not, follow the instructions here: [https://julialang.org/downloads/platform/](https://julialang.org/downloads/platform/).

## Install PackageToBeShipped
### Activate an environment. 
For beginners of Julia, it might be a good idea to 

Press ] to activate the package REPL. 
type 

```julia
activate EnvironmentName
```

where EnvironmentName is the name you want to give to your environment. 
Next install PackageToBeShipped to your environment: 

```julia
add PackageToBeShipped
```

If you want you can install more packages to your environment. For instance, most likely, you will need, CSV, DataFrames,
GZip, CodecZlib, StatsBase, LinearAlgebra, Distributions, Distributed

```julia
] add CSV DataFrames GZip CodecZlib StatsBase LinearAlgebra Distributions Distributed
```

Read more about Julia environments here [https://pkgdocs.julialang.org/v1/environments/](https://pkgdocs.julialang.org/v1/environments/).

## Data format. 

Our program can work with genetic data in the DICE-2 format. See [https://github.com/grenaud/dice](https://github.com/grenaud/dice?tab=readme-ov-file#2-pop-method-input-data-format). 

This is a CSV-file with four columns, in this order: number of ancestral reads, number of derived reads, the frequency of the anchor population, and in how many loci this combination of reads and frequency occurs. So the first two columns contain non-negative integers, the third column a real number between 0.0 and 1.0 and the last a positive integer. 



## Using the sampler. 

You can use the example script in [example.org](example.org).

This are the steps that you need to take. 
1. Put your data in the DICE-2 format, and save it as a (zipped) CSV file. 
2. Specify a prior on n, (τC, τA) and ϵ. For example you could use a uniform prior for all three:
```julia
using Distributions # Defines probabililty distributions. 
prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]
prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}
prioronϵ = Uniform(0, 0.5) # uniform prior on interval[0,0.5]
```
3. Specify how many mcmc samples you want and how many chains you want to sample. For example, if you want 4 chains with each 100'000 samples you can specify 
```julia
nchains = 4 
nsteps = 100_000
```
4. Next you need to set the path to your DICE file. For example 
```
dicefile = "path/to/dicefile.csv" # UNIX-systems
dicefile = "path\\to\\dicefile.csv" # Windows-sytems. 
```
where "path/to/dicefile.csv" is the relative path to dicefile.csv from your working directory. You can find the working directory as follows: Press ";" to come in the shell REPL.  
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
using PackageToBeShipped, Distributions, CSV  
using 




prioronτCτA = Product([Uniform(), Uniform()]) # uniform prior on [0, 1] x [0, 1]

prioronn = DiscreteUniform(1, 10) # uniform prior on {1, 2, ..., 10}

prioronϵ = Uniform(0, 0.5) # uniform prior on [0, 0.5]

nchains = 4 

nsteps = 100_000

dicefile = "path/to/dicefile.csv"

resultscsvfile = "path/to/results.csv"

resultsdataframe = MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, dicefile)

CSV.write(resultscsvfile, resultsdataframe)
```

Congratulations! You got your mcmc samples. They are stored as a DataFrame and on the disk as a CSV file. You can use the CSV file for further analysis. 
