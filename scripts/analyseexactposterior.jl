# Created 13 - 10 - 2023
# Set of Julia instructions to analyse the exact posterior.
# Meant to use interactively. 

using CSV, DataFrames, Plots 

ancientpops = ["CHBS", "YRIS", "CEUS", "KARS"]
anchorpops = ["YRI", "CEU", "CHB", "KAR"]
truenumindvec = [1, 2, 4, 8]
allcombinations = [(anchorpop, ancientpop, truenumind) for anchorpop in anchorpops for ancientpop in ancientpops for truenumind in truenumindvec] # 64 elementen. 

trueϵ = 0.001 # true error rate used in the simulations. 
coverage = 64 # the coverage is 64 everywhere. 
τrange = 0.001:0.001:1.0

allcombinations2 = allcombinations[[38, 63, 36, 57]]
allcombinations3 = setdiff(allcombinations, allcombinations2)
for (anchorpop, ancientpop, truenumind) in allcombinations3
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "exactposterior", string("Madsmodel_ancient_pop_exact_posterior_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    df = CSV.read(zippedcsvfile, DataFrame)
    plotfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "plots", "Madsmodel", string("Madsmodel_plot_exact_posterior_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".png"))
    p = surface(df.τCvalue, df.τAvalue, df.logjointprob)
    savefig(p, plotfile)
end 


#=
using Plots
gr()
p = surface(df.τCvalue[begin:1000:end], df.τAvalue[begin:1000:end], df.logjointprob[begin:1000:end]);
savefig()
=#

#=
f(x) = in(x, 0.00:0.01:1)

df2 = df[f.(df.τCvalue) .& f.(df.τAvalue), :]
using Plots
gr()
p = surface(df2.τCvalue, df2.τAvalue, df2.logjointprob)
=#