csvversion, distributionsversion, statsbaseversion = ARGS  

f = open("compatibleversions.txt", "a")

using Pkg 
Pkg.activate("Testversions")
try
    Pkg.add(name="CSV", version=csvversion)
    Pkg.add(name="DataFrames", version="1")
    Pkg.add(name="Distributions", version=distributionsversion)
    Pkg.add(name="LinearAlgebra", version="1")
    Pkg.add(name="StaticArrays", version="1")
    Pkg.add(name="StatsBase", version=statsbaseversion)
    println(f, "CSV $csvversion Distributions $distributionsversion StatsBase $statsbaseversion succeeded")
catch 
    println(f, "CSV $csvversion Distributions $distributionsversion StatsBase $statsbaseversion did not succeed")
end 

close(f)

Pkg.status()
println("finished")