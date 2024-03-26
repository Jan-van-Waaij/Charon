# analyseposterior.jl
# Author: Jan van Waaij
# Created: 6 March 2024
# Last update: 6 March 2024
#
# Purpose:
# This is a collection of Julia commands to analyse the MCMC samples.
#
# Installation:
# You need the GZip, DataFrames, Plots and CSV packages. 
# Please read the Julia package install instructions in installinstructions.txt.
# 
# Usage:
# open a Julia terminal and use snippets of the code interactively. 
#
# Requirement:
# a csv.gz file with MCMC samples from /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/mcmcsamples

# UNKNOWN N and with errors:

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      CSV, # for CSV files.
      StatsBase, # for the mean function
      Plots, 
      StatsPlots,
      MCMCChains

ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec] 
nsteps = 400_000 # number of steps in mcmc chain 
nchains = 4 # number of chains. 
trueϵ = 0.001 # true error rate used in the simulations. 
coverage = 64 # the coverage is 64 everywhere. 

anchorpop = "CHB"
truenumind = 2
ancientpop = "YRIS"

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
zippedcsvfile = csvfile * ".gz"

zippedcsvfile = "data/mcmcsamples/Madsmodel_ancient_pop_KARS_anchor_pop_KAR_true_num_ind_4_factor_1.0.csv.gz"

df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])

gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])

# make plots
#allcombinations = allcombinations[[38, 63, 36, 57]]

stdsparameters = [(10.0^(-3), 0.001), (10.0^(-4.5), 0.0001)]
allcombinations2 = allcombinations[[38, 63, 36, 57]]
allcombinations3 = [(tuple1..., tuple2...) for tuple1 in allcombinations2 for tuple2 in stdsparameters] 
for (anchorpop, ancientpop, truenumind, stdτ, stdϵ) in allcombinations3
    if (stdτ, stdϵ) != stdsparameters[2]
        continue 
    end 
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_Special_MCMC_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, "_stdtau_", stdτ, "_stdeps_", stdϵ, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    # chainarray = cat((Matrix(gdf[k][nsteps÷2:nsteps, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    # chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
    chainarray = cat((Matrix(gdf[k][nsteps÷2:nsteps, [:τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    chain = Chains(chainarray, [:τCsample, :τAsample, :ϵsample])
    p=plot(chain)
    display(p)
end 


allcombinations2 = allcombinations[[38, 63, 36, 57]]
temperatures = [200, 400, 1000, 10_000.0]
allcombinations3 = [(tuple..., temperature) for tuple in allcombinations2 for temperature in temperatures]
for (anchorpop, ancientpop, truenumind, temperature) in allcombinations3
    if temperature != 10_000.0 
        continue 
    end 
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, "_temperature_", temperature, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    # chainarray = cat((Matrix(gdf[k][nsteps÷2:nsteps, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    # chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
    chainarray = cat((Matrix(gdf[k][nsteps÷2:nsteps, [:τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    chain = Chains(chainarray, [:τCsample, :τAsample, :ϵsample])
    p=plot(chain)
    display(p)
end 

# calculate acceptancerate
for (anchorpop, ancientpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    for chain in 1:4
        rate = mean(gdf[chain].:accepted[end÷8:end÷4])
        println("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " rate ", rate)
    end 
end 

# calculate acceptancerate
stdsparameters = [(10.0^(-3), 0.001), (10.0^(-4.5), 0.0001)]
allcombinations2 = allcombinations[[38, 63, 36, 57]]
allcombinations3 = [(tuple1..., tuple2...) for tuple1 in allcombinations2 for tuple2 in stdsparameters]
for (anchorpop, ancientpop, truenumind, stdτ, stdϵ) in allcombinations3
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_Special_MCMC_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, "_stdtau_", stdτ, "_stdeps_", stdϵ, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    for chain in 1:4
        rate = mean(gdf[chain].:accepted[end÷2:end])
        println("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " std τ ", stdτ, " std ϵ ", stdϵ,  " rate ", rate)
    end 
end 

# calculate acceptancerate
allcombinations2 = allcombinations[[38, 63, 36, 57]]
temperatures = [200, 400, 1000, 10_000.0]
allcombinations3 = [(tuple..., temperature) for tuple in allcombinations2 for temperature in temperatures]
for (anchorpop, ancientpop, truenumind, temperature) in allcombinations3
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, "_temperature_", temperature, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    for chain in 1:4
        rate = mean(gdf[chain].:accepted[end÷2:end])
        println("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " temperature ", temperature, " rate ", rate)
    end 
end 

isconstant(x) = length(unique(x)) == 1 
# constant n for all chains? 
constantchainsdisagree = Bool[]
for (anchorpop, ancientpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    chainsagree = isconstant(vcat([gdf[k].:nsample[end÷2:end] for k in eachindex(gdf)]...))
    everyseperatechainisconstant = all([isconstant(gdf[k].:nsample[end÷2:end]) for k in eachindex(gdf)])
    push!(constantchainsdisagree, everyseperatechainisconstant & !chainsagree)
    #println("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " chains agree? ", chainsagree, " every chain is constant? ", everyseperatechainisconstant) 
end 

constantchainsdisagree = Bool[]
chainsagreevec = Bool[]
for (anchorpop, ancientpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    chainsagree = isconstant(vcat([gdf[k].:nsample[end÷2:end] for k in eachindex(gdf)]...))
    push!(chainsagreevec, chainsagree)
    everyseperatechainisconstant = all([isconstant(gdf[k].:nsample[end÷2:end]) for k in eachindex(gdf)])
    push!(constantchainsdisagree, everyseperatechainisconstant & !chainsagree)
    #println("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " chains agree? ", chainsagree, " every chain is constant? ", everyseperatechainisconstant) 
end 

f(x) = mean(x[end÷2:100:end])
g(x) = median(x[end÷2:100:end])

minn(x) = minimum(x[end÷2:100:end])
maxn(x) = maximum(x[end÷2:100:end])

lowerbound(x) = quantile(x[end÷2:100:end], 0.025)
upperbound(x) = quantile(x[end÷2:100:end], 0.975)

#cdf= combine(gdf, :nsample => f => mean)

grid = 1:10_000:80_001  # vcat(1:80_000:400_000, 400_000)
ranges = [grid[n]:grid[n+1] for n in eachindex(grid)[1:end-1]]

acceptance(x, range) = mean(x[range])

cdf = combine(gdf, :accepted .=> [x -> acceptance(x, range) for range in ranges] .=> [Symbol(:accrange, n) for n in 1:length(ranges)])

allcdf2 = DataFrame()
for (anchorpop, ancientpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    cdf = combine(gdf, :accepted .=> [x -> acceptance(x, range) for range in ranges] .=> [Symbol(:accrange, n) for n in 1:length(ranges)])
    append!(allcdf2, cdf)
end 


cdf= combine(gdf, :nsample => lowerbound => Symbol("n 0.025%"),
                     :nsample => f => :nmean,
                     :nsample => upperbound => Symbol("n 0.975%"),
                     :τCsample => lowerbound => Symbol("τC 0.025%"),
                     :τCsample => f => :meanτC,
                     :τCsample => g => :medianτC,
                     :τCsample => upperbound => Symbol("τC 0.975%"),
                     :τAsample => lowerbound => Symbol("τA 0.025%"),
                     :τAsample => f => :meanτA,
                     :τAsample => g => :medianτA,
                     :τAsample => upperbound => Symbol("τA 0.975%"),
                     :ϵsample => lowerbound => Symbol("ϵ 0.025%"),
                     :ϵsample => f => :meanϵ,
                     :ϵsample => upperbound => Symbol("ϵ 0.975%"),
                     :accepted => f => :acceptanceratelasthalf)

#show(cdf[cdf.coverage .== 7, [ :coverage, :anchorpop, :truenumind, :nmean, :trueτC, :meanτC, :trueτA, :meanτA, :trueϵ, :meanϵ]], allrows=true)

show(cdf[:, [ :anchorpopulation, :truenumind, :nmean, :trueτC, :medianτC, :trueτA, :medianτA, :trueϵ, :meanϵ]], allrows=true)

show(cdf[:, [ :anchorpopulation, :truenumind, :nmean, :trueτC, :meanτC, Symbol("τA 0.025%"), :trueτA, Symbol("τA 0.975%"), :meanτA, :trueϵ, :meanϵ]], allrows=true)

pdf= cdf[:, [:coverage, :anchorpop, :truenumind, Symbol("n 0.025%"), Symbol("n 0.975%"), :trueτC, Symbol("τC 0.025%"), Symbol("τC 0.975%"), :trueτA, Symbol("τA 0.025%"), Symbol("τA 0.975%"), Symbol("ϵ 0.025%"), Symbol("ϵ 0.975%")]]

show(pdf[pdf.coverage .== 7, :], allrows=true)

show(cdf[:, [:acceptanceratelasthalf, :truenumind, :nmean, :trueτA, :meanτA]], allrows=true)
        

show(cdf[:, [:acceptanceratelasthalf, #:temperature, 
:truenumind, :nmean]], allrows=true)

show(cdf[:, [:acceptanceratelasthalf, #:temperature, 
:truenumind, Symbol("n 0.025%"), Symbol("n 0.975%")]], allrows=true)


confint = cdf[:, [Symbol("τA 0.025%"), :trueτA , Symbol("τA 0.975%"),]]

confint.contained = confint[:, 1] .≤ confint[:, 2] .≤ confint[:, 3]

show(confint, allrows=true)




show(cdf[:, [:acceptanceratelasthalf, #:temperature,
 :trueτC, :meanτC]], allrows=true)

show(cdf[:, [:acceptanceratelasthalf, #:temperature,
 :trueτA, :meanτA]], allrows=true)


show(cdf[:, [:acceptanceratelasthalf, #:temperature, 
:trueϵ, :meanϵ]], allrows=true)

pdf = cidf[:, ["NeAnchorPop", "NeAncientPop", "jointtimeyears",  "n 0.025%", "truenumind", "n 0.975%",   "τC 0.025%", "trueτC", "τC 0.975%",  "τA 0.025%", "trueτA", "τA 0.975%", "acceptanceratelasthalf"]]

pdf = cdf[:, [:acceptanceratelasthalf, :truenumind, :nmean, :trueτC, :meanτC, :trueτA, :meanτA, :trueϵ, :meanϵ]]

pdf = cdf[:, [:acceptanceratelasthalf, :truenumind, :nmean, :trueτC, :meanτC, :trueτA, :meanτA, :trueϵ, :meanϵ]]


a(x::Int) = x
a(x::Float64) = round(x, sigdigits=3)

pdf2 = a.(pdf)

using Latexify # make nice latex code from 
# julia objects

latexify(pdf2, env = :table)



sacdf2 = acdf2[:, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :chainid, :truenumind, :minn, :meann, :maxn, :stdn, :trueτC, :meanτC, :stdτC, :trueτA, :meanτA, :stdτA, :acceptanceratelasthalf]]

a(x::Int) = x
a(x::Float64) = round(x, digits=2)

sacdf3 = a.(sacdf2)

show(sacdf3, allrows=true)


latexify(sacdf3, env = :table)


using Plots,
      StatsPlots,
      MCMCChains

nsteps = 100_000
nchains = 4
#picindex = 0
for numdataset in 1:length(gdf)÷nchains
    gdfindex = (numdataset-1)*nchains # k = 1,..., nchains 
    
    chainarray = cat((Matrix(gdf[gdfindex+k][1:nsteps, 
    [:nsample, 
    :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    #global picindex += 1
    chain = Chains(chainarray, [:nsample, 
    :τCsample, :τAsample, :ϵsample])
    display(plot(chain, #title = string(gdf[gdfindex+1].temperature[1])
    ))
    #sleep(4)
end 
println("hi")


#=
Conclusion:
With unknown N, and errors, starting in a starting value choosen by the prior, can have really bad mixing and non-convergence properties. However, starting in the true value leads to perfect mixing. 
=# 

# UNKNOWN N:

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      Plots, # for plotting,
      StatsPlots, # for displaying MCMC chains. 
      CSV, # for CSV files.
      StatsBase, # for the mean function
      MCMCChains, # for analysing MCMC chains on things like convergence.
      Latexify # make nice latex code from 
      # julia objects
      
csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "mcmcwithlocalmaximacoverage5.csv.gz")
# data/mcmcsamples/allsamplesunknownnwithtrueinitialvalue.csv.gz
#names = [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :truenumind, :τCsample, :τAsample, :accepted, :chainid]

df = CSV.read(csvfile, DataFrame) #, header = names)

df2 = df[df.coverage .== 100 , :]

gdf = groupby(df, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :truenumind, :coverage, :localmaxnumber, :chainid])

df3 = vcat((gdf[k][end÷2:100:end,:] for k in 1:length(gdf))...)

gdf2 = groupby(df3, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :truenumind])


f(x) = mean(x[end÷2:100:end])

g(x) = mean(x[end÷2:100:end])

h(x) = std(x[end÷2:100:end])

minn(x) = minimum(x[end÷2:100:end])
maxn(x) = maximum(x[end÷2:100:end])

lowerbound(x) = quantile(x[end÷2:100:end], 0.025)
upperbound(x) = quantile(x[end÷2:100:end], 0.975)


cdf= combine(gdf, :nsample => lowerbound => Symbol("n 0.025%"),
                     :nsample => mean => :nmean,
                     :nsample => upperbound => Symbol("n 0.975%"),
                     :τCsample => lowerbound => Symbol("τC 0.025%"),
                     :τCsample => mean => :meanτC,
                     :τCsample => upperbound => Symbol("τC 0.975%"),
                     :τAsample => lowerbound => Symbol("τA 0.025%"),
                     :τAsample => mean => :meanτA,
                     :τAsample => upperbound => Symbol("τA 0.975%"),
                     :accepted => f => :acceptanceratelasthalf)

pdf = cidf[:, ["NeAnchorPop", "NeAncientPop", "jointtimeyears",  "n 0.025%", "truenumind", "n 0.975%",   "τC 0.025%", "trueτC", "τC 0.975%",  "τA 0.025%", "trueτA", "τA 0.975%", "acceptanceratelasthalf"]]

a(x::Int) = x
a(x::Float64) = round(x, digits=2)

pdf2 = a.(pdf)

latexify(pdf2, env = :table)



sacdf2 = acdf2[:, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :chainid, :truenumind, :minn, :meann, :maxn, :stdn, :trueτC, :meanτC, :stdτC, :trueτA, :meanτA, :stdτA, :acceptanceratelasthalf]]

a(x::Int) = x
a(x::Float64) = round(x, digits=2)

sacdf3 = a.(sacdf2)

show(sacdf3, allrows=true)


latexify(sacdf3, env = :table)

numdataset=1
nchains = 4
picindex = 0
for numdataset in 1:length(gdf)÷nchains
    gdfindex = (numdataset-1)*nchains # k = 1,..., nchains 
    
    chainarray = cat((Matrix(gdf[gdfindex+k][25_001:100_000, #10_000:24_999, #25_000:100_000, 
    [:nsample, :τCsample, :τAsample]]) for  k in 1:nchains)..., dims=3)
    global picindex += 1
    chain = Chains(chainarray, [:nsample, :τCsample, :τAsample])
    display(plot(chain))
    sleep(4)
end 


numdataset=1
nchains = 4
picindex = 0
for numdataset in 1:length(gdf)÷nchains
    gdfindex = (numdataset-1)*nchains # k = 1,..., nchains 
    
    chainarray = cat((Matrix(gdf[gdfindex+k][75_001:100_000, #10_000:24_999, #25_000:100_000, 
    [:nsample, :τCsample, :τAsample]]) for  k in 1:nchains)..., dims=3)
    global picindex += 1
    info = (picindex = picindex, NeAnchorPop = gdf[gdfindex+1][1, :NeAnchorPop], NeAncientPop = gdf[gdfindex+1][1, :NeAncientPop], jointtimeyears = gdf[gdfindex+1][1, :jointtimeyears], truenumind = gdf[gdfindex+1][1, :truenumind], acceptancerate = sacdf3[gdfindex+1:gdfindex+nchains, :acceptanceratelasthalf])
    chain = Chains(chainarray, [:nsample, :τCsample, :τAsample], info = info)
    @show info 
    display(plot(chain))
    sleep(4)
end 

function displaychains(numdataset)
    gdfindex = (numdataset-1)*nchains # k = 1,..., nchains 
    
    chainarray = cat((Matrix(gdf[gdfindex+k][75_001:100_000, #10_000:24_999, #25_000:100_000, 
    [:nsample, :τCsample, :τAsample]]) for  k in 1:nchains)..., dims=3)
    global picindex += 1
    info = (picindex = picindex, NeAnchorPop = gdf[gdfindex+1][1, :NeAnchorPop], NeAncientPop = gdf[gdfindex+1][1, :NeAncientPop], jointtimeyears = gdf[gdfindex+1][1, :jointtimeyears], trueτC = round(gdf[gdfindex+1][1, :trueτC], digits=2), trueτA = round(gdf[gdfindex+1][1, :trueτA], digits = 2), truenumind = gdf[gdfindex+1][1, :truenumind], acceptancerate = sacdf3[gdfindex+1:gdfindex+nchains, :acceptanceratelasthalf])
    chain = Chains(chainarray, [:nsample, :τCsample, :τAsample], info = info)
    @show info 
    display(plot(chain))
end 

println("ready")

# KNOWN N:

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      Plots, # for plotting,
      StatsPlots, # for displaying MCMC chains. 
      CSV, # for CSV files.
      StatsBase, # for the mean function
      MCMCChains, # for analysing MCMC chains on things like convergence.
      Latexify # make nice latex code from 
      # julia objects

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", "allsamplesfixedndifferentcoveragesanddifferentnumindividualsmanyindividuals.csv.gz")

#names = [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :truenumind, :τCsample, :τAsample, :accepted, :chainid]

df = CSV.read(csvfile, DataFrame) #, header = names)

gdf = groupby(df, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :truenumind, :coverage, :chainid])

f(x) = mean(x[end÷2:end])

g(x) = mean(x[end÷2:100:end])

h(x) = std(x[end÷2:100:end])

acdf2 = combine(gdf, :τCsample => g => :meanτC, :τCsample => h => :stdτC, :τAsample => g => :meanτA, :τAsample => h => :stdτA, :accepted => f => :acceptanceratelasthalf)

sacdf2 = acdf2[:, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :chainid, :coverage, :truenumind,  :trueτC, :meanτC, :stdτC, :trueτA, :meanτA, :stdτA, :acceptanceratelasthalf]]

a(x::Int) = x
a(x::Float64) = round(x, digits=2)

sacdf3 = a.(sacdf2)

show(sacdf3, allrows=true)


latexify(sacdf3, env = :table)




cdf = combine(gdf, [:loglik, :numind, :τC, :τA] 
                    => ((x,y...) -> reshape([z[argmax(x)] for z in y], 1, :)) 
                    => [:maxliknumind, :maxlikτC, :maxlikτA])

# shuffle the columns. 
scdf = cdf[:, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :truenumind, :maxliknumind, :trueτC, :maxlikτC, :trueτA, :maxlikτA]]

# show all columns. 
show(scdf, allrows=true, eltypes=false)


round.(scdf, digits=2)
# analyse the chains with MCMCChains package. 

numdataset=1 
nchains = 4
picindex = 0
for numdataset in 1:length(gdf)÷nchains
    gdfindex = (numdataset-1)*nchains # k = 1,..., nchains 
    
    chainarray = cat((Matrix(gdf[gdfindex+k][1:2_000, #25_000:100_000, 
    [:τCsample, :τAsample]]) for  k in 1:nchains)..., dims=3)
    global picindex += 1
    info = (picindex = picindex, NeAnchorPop = gdf[gdfindex+1][1, :NeAnchorPop], NeAncientPop = gdf[gdfindex+1][1, :NeAncientPop], jointtimeyears = gdf[gdfindex+1][1, :jointtimeyears], truenumind = gdf[gdfindex+1][1, :truenumind], coverage = gdf[gdfindex+1][1, :coverage], acceptancerate = sacdf3[gdfindex+1:gdfindex+nchains, :acceptanceratelasthalf])
    chain = Chains(chainarray, [:τCsample, :τAsample], info = info);
    @show info 
    display(plot(chain))
    sleep(4)
end 
println("ready")
