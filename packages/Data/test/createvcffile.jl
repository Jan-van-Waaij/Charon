using StatsBase, GZip

AA = "1|1"
Aa = "1|0"
aA = "0|1"
aa = "0|0"
fiftypercentderived = vcat([repeat([x], 100) for x in [AA, Aa, aA, aa]]...)
len = length(fiftypercentderived)

numlines = 10_000

file = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "packages", "Data", "test", "fiftypercentderived.vcf.gz")
openedfile = GZip.open(file, "w")

for _ in Base.OneTo(numlines)
    genotypesvec = sample(fiftypercentderived, len, replace = false)
    genotypesvecwithtabs = vcat([[a, "\t"] for a in genotypesvec]...)[begin:end-1]
    genotypesstring = string(genotypesvecwithtabs...) 
    println(openedfile, genotypesstring)
end 

close(openedfile)