csvversion, distributionsversion, statsbaseversion = ARGS  

f = open("compatibleversionsonlynumbers.txt", "a")

using Pkg 
Pkg.activate("Testversionsonlynumbers")
try
    Pkg.add(name="CSV", version=csvversion)
    Pkg.add(name="DataFrames", version="1")
    Pkg.add(name="Distributions", version=distributionsversion)
    Pkg.add(name="LinearAlgebra", version="1")
    Pkg.add(name="StaticArrays", version="1")
    Pkg.add(name="StatsBase", version=statsbaseversion)
    println(f, "$csvversion $distributionsversion $statsbaseversion")
catch 
end 

close(f)

Pkg.status()
println("finished")