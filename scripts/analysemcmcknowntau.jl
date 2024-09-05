# analyseposterior.jl
# Author: Jan van Waaij
# Created: 9 May 2024
# Last update: 9 May 2024
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




forestplot([[0.2,1.2], [0.8,1.4], [0.6, 0.8], [1.2, 1.6]], 
       metric = [1.0, 1.2, 0.7, 1.3], source = ["A", "B", "C", "D"],
       summary = Dict(:ci =>[0.8, 1.1], :est => 0.95, :markershape => :rtriangle), logscale = false, refline=0.25, sourcelabel="coverage", metriclabel="estimate")

using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      CSV, # for CSV files.
      StatsBase, # for the mean function
      Plots, # plotting. 
      StatsPlots, 
      MCMCChains, # plotting of mcmc chains. 
      Distributed, # distributed programming
      Latexify, # make latex code from Julia objects. 
      ForestPlot

forestplot = ForestPlot.forestplot

dir = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/images/empirical_data/"
exename = `nice -19 /home/ctools/julia-1.6.7/bin/julia`
exeflags = ["--project=$dir", "--optimize=3"]
      
nodes = [
   #= ("compute01", 30), # 50 percent used. MAX CAP FOR ME: 31 ASO MODE 32 
    ("node05", 6), #6), # 2 cores in use. MAX CAP FOR ME:5 ASO MODE 6
    ("node03", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node04", 6), # not used MAX CAP FOR ME: 5 ASO MODE 6
    ("node12", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4 
    ("node01", 6), # not used. MAX CAP FOR ME: 5 ASO MODE 6
    ("node02", 4), # not used MAX CAP FOR ME: 3 ASO MODE 4
    ("node06", 6), # 2 core in use. MAX CAP FOR ME: 5 ASO MODE 6
    ("node07", 13), #14), # 16), # fully used MAX CAP FOR ME: 13 ASO MODE 16
    ("node13", 4), # 2 used. MAX CAP FOR ME: 3 ASO MODE 4
    ("node11", 10), #10), # fully used MAX CAP FOR ME: 8 ASO MODE 10
  =#
    ("node14", 16) #, # not used MAX CAP FOR ME: 3 ASO MODE 4
#    ("node08", 5), # not used MAX CAP FOR ME: 4 ASO MODE 5
]
sum(node[2] for node in nodes)


addprocs(nodes, dir=dir, exename=exename, exeflags=exeflags)

@everywhere ENV["GKSwstype"]="100"
@everywhere using GZip, # to work with zipped files,
DataFrames, # to work with DataFrames,
CSV, # for CSV files.
StatsBase, # for the mean function
Plots, # plotting. 
StatsPlots, 
MCMCChains, # plotting of mcmc chains. 
Distributed # distributed programming


numberofindividualsinthesedimentvec = 4:4
coverages = 15:-1:10
allcombinations4 = [(numind, cov) for numind in numberofindividualsinthesedimentvec for cov in coverages]

numberofindividualsinthesedimentvec = 6:6
coverages = vcat(100, 50, 40, 30:-1:20, 15, 10)
allcombinations = [(numind, cov) for numind in numberofindividualsinthesedimentvec for cov in coverages]


numberofindividualsinthesedimentvec = 7:7
coverages = [100, 50:-1:30..., 20, 15, 10]
allcombinations = [(numind, cov) for numind in numberofindividualsinthesedimentvec for cov in coverages]



numberofindividualsinthesedimentvec = 8:8
coverages = vcat(100, 80:-1:70, 60:-10:20, 15, 10)
allcombinations = [(numind, cov) for numind in numberofindividualsinthesedimentvec for cov in coverages]

allcombinations = vcat(allcombinations4, allcombinations6, allcombinations7, allcombinations8)


numberofindividualsinthesedimentvec = 4:4
coverages = 15:-1:10
allcombinations = [(numind, cov) for numind in numberofindividualsinthesedimentvec for cov in coverages]


nsteps = 100_000 # number of steps in mcmc chain 
nchains = 4 # number of chains. 
trueϵ = 0.003 # true error rate used in the simulations. 

results = DataFrame(
    truenumind = Int[],
    coverage = Int[], 
    n_low = Float64[],
    n_median = Float64[],
    n_high = Float64[],
    τC_low = Float64[],
    τC_median = Float64[],
    τC_high = Float64[],
    τA_low = Float64[],
    τA_median = Float64[],
    τA_high = Float64[],
    ϵ_low = Float64[],
    ϵ_median = Float64[],
    ϵ_high = Float64[],
)

allind = 1:8
allcov = 1:100
allcombinations = [(numind, cov) for numind in allind for cov in allcov]

for (numberofindividualsinthesediment, coverage) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", coverage, "_MCMC.csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        println(zippedcsvfile, " does not exist.")
        continue
    end 
    # the following code checks whether csvfile is not yet added to results.
    if !isempty(filter([:truenumind, :coverage] => (x,y) -> x == numberofindividualsinthesediment && y == coverage , results))
        continue 
    end 
    df = CSV.read(zippedcsvfile, DataFrame, header = ["trueϵ", "truenumind", "coverage", "nsample", "τCsample", "τAsample", "ϵsample", "accepted", "logjointprob", "chainid"])
    # sometimes something goes wrong with the header, then we try it again.
    # but this time without adding the header names.  
    if eltype(df.trueϵ) <: AbstractString
        df = CSV.read(zippedcsvfile, DataFrame)
    end 
    # remove the first half and thin the sample. 
    df = df[end÷2:10:end,:]
    # separate the dataframe by chain. 
    gdf = groupby(df, [:truenumind, :coverage, :trueϵ, :chainid])
    # calculate the credible intervals.
    n_low, n_median, n_high = quantile(df.nsample, [0.025, 0.5, 0.975])
    τC_low, τC_median, τC_high = quantile(df.τCsample, [0.025, 0.5, 0.975])
    τA_low, τA_median, τA_high = quantile(df.τAsample, [0.025, 0.5, 0.975])
    ϵ_low, ϵ_median, ϵ_high = quantile(df.ϵsample, [0.025, 0.5, 0.975])
    # add the results to the dataframe. 
    push!(results, (truenumind = numberofindividualsinthesediment, 
                    coverage = coverage, 
                    n_low = n_low, 
                    n_median = n_median, 
                    n_high = n_high, 
                    τC_low = τC_low, 
                    τC_median = τC_median, 
                    τC_high = τC_high, 
                    τA_low = τA_low, 
                    τA_median = τA_median, 
                    τA_high = τA_high, 
                    ϵ_low = ϵ_low, 
                    ϵ_median = ϵ_median, 
                    ϵ_high = ϵ_high))
    println("Finished with ", csvfile)
    println(results[end, :])
end 

results_sorted = sort(results, [:truenumind, :coverage])

for numind in 1:8 
    resultsfixedn = filter(:truenumind => x -> x .== numind, results)

    crediblesets = [[x, y] for (x,y) in zip(resultsfixedn.:τC_low, resultsfixedn.:τC_high)]

    source = resultsfixedn.:coverage

    metric = resultsfixedn.:τC_median 

    myforestplot(crediblesets, source = source, metric = metric, logscale = false, refline=0.25, sourcelabel="coverage", metriclabel="estimate")

    file = "data/images/subsampling/credibleintervals/subsampling_number_of_individuals=$(numind)_credibleinterval_tauC.png"

    savefig(file)

    crediblesets = [[x, y] for (x,y) in zip(resultsfixedn.:τA_low, resultsfixedn.:τA_high)]

    source = resultsfixedn.:coverage

    metric = resultsfixedn.:τA_median 

    myforestplot(crediblesets, source = source, metric = metric, logscale = false, refline=0.25, sourcelabel="coverage", metriclabel="estimate")

    file = "data/images/subsampling/credibleintervals/subsampling_number_of_individuals=$(numind)_credibleinterval_tauA.png"

    savefig(file)
end 

for numind in 1:8 
    resultsfixedn = filter(:truenumind => x -> x .== numind, results)

    crediblesets = [[x, y] for (x,y) in zip(resultsfixedn.:n_low, resultsfixedn.:n_high)]

    source = resultsfixedn.:coverage

    metric = resultsfixedn.:n_median 

    myforestplot(crediblesets, source = source, metric = metric, logscale = false, refline=0.25, sourcelabel="coverage", metriclabel="estimate")

    file = "data/images/subsampling/credibleintervals/subsampling_number_of_individuals=$(numind)_credibleinterval_n.png"

    savefig(file)
end 

function test(ci)
    maxx = round((maximum(x -> x[2], ci) - 0.1)/2, digits = 1)*2 + 0.2
    minx = round((minimum(x -> x[1], ci) - 0.1)/2, digits = 1)*2 
    xlims = (minx, maxx)
    ticks = collect(minx:0.2*floor((maxx - minx)/0.2/4):maxx)
    return length(ticks)
end


function myforestplot(ci; sourcelabel = "Source:", metriclabel = "OR", cilabel = "CI95%", 
    source = nothing, metric = nothing, printci = false, refline::Union{Real, Bool} = 1,
    summary = nothing, logscale = true, cimsz = -1, cimszwts = nothing, size = (600, 400), ps = 10, kwargs...)

    #size=(800,400)
    lines = length(ci) + 3
    if lines <= 6  lines = 8 end
    ylims = (0, lines)
    sty   = lines - 1
    

    if logscale 
        func = exp
    else
        func = identity
    end

    maxx = round((maximum(x -> x[2], ci) - 0.1)/2, digits = 1)*2 + 0.2
    minx = round((minimum(x -> x[1], ci) - 0.1)/2, digits = 1)*2 
    xlims = func.((minx, maxx))
    ticks = collect(range(minx, maxx, length=5))
    
    p = plot(size = size)

    for i = 1:length(ci)
        plot!(p, func.(ci[i]), fill(sty - i, 2), markershape = :vline, linecolor = :blue, markercolor = :blue)
    end

    if cimsz < 0 && (!isnothing(metric) || !isnothing(cimszwts)) 
        if isnothing(cimszwts) cimszwts = metric end
        minm = minimum(cimszwts)
        c = (maximum(cimszwts) - minm)/3 
        cimsz = @. (cimszwts - minm) / c + 3 
    end

    if !isnothing(metric)
        my = sty .- collect(1:length(metric))
        plot!(func.(metric), my, seriestype= :scatter, markershape = :rect, markercolor = :gray, markeralpha = 0.8, markersize = cimsz, legend = false)
    end

    # Print summary CI
    if !isnothing(summary) && haskey(summary, :ci) && haskey(summary, :est)
        if haskey(summary, :markershape)
            ms = summary[:markershape]
        else
            ms = :diamond
        end
        if haskey(summary, :markersize)
            msz = summary[:markersize]
        else
            msz = 5
        end
        plot!(p, func.(summary[:ci]), fill(1, 2), markershape = :vline, linecolor = :red, markercolor = :red)
        plot!([func(summary[:est])], [1], seriestype= :scatter, markershape = ms, markercolor = :red, legend = false, markersize = msz)
    end
    if !isnothing(summary) && haskey(summary, :vline) && !(summary[:vline] === false) # Print Vlines
        if summary[:vline] === true && haskey(summary, :est) 
            plot!(p, func.([summary[:est], summary[:est]]),[0, sty], line = (:dash, 1.5), linealpha = 0.7, linecolor = :red)
        else 
            for v in summary[:vline]
                plot!(p, func.([v, v]),[0, sty], line = (:dash, 1.5), linealpha = 0.7, linecolor = :red)
            end
        end 
    end
    if !(refline === false)
        plot!(p, func.([refline, refline]),[0, sty]; linetype=:line, linealpha = 0.7,  linecolor = :black, line = (:dot, 1.5))
    end
    plot!(p; ylims = ylims, xlims = xlims, yshowaxis = false, ticks = (func.(ticks), round.(ticks, sigdigits=3)), kwargs...) # (@. round(func(ticks), digits=2), ticks))

    if !isnothing(source)
        tp = plot(showaxis = false, xlims=(0, 2.5), ylims = ylims, size = size)
        t = Plots.text(sourcelabel, halign = :left, family = "Palatino Bold", pointsize = ps)
        annotate!(tp, 0, sty, t)
        for i = 1:length(source)
            t = Plots.text(source[i], halign = :left, family = "Palatino", pointsize = ps)
            annotate!(tp, 0, sty - i, t)
        end

        if !isnothing(metric) 
            t = Plots.text(metriclabel, halign = :left, family = "Palatino Bold", pointsize = ps)
            annotate!(tp, 1.6, sty, t)
            for i = 1:length(metric)
                t = Plots.text(string(round(metric[i], digits = 3)), halign = :left, family = "Palatino", pointsize = ps)
                annotate!(tp, 1.8, sty - i, t)
            end
        end

        if printci
            t = Plots.text(cilabel, halign = :left, family = "Palatino Bold", pointsize = ps)
            annotate!(tp, 2.3, sty, t)
            for i = 1:length(ci)
                t = Plots.text(string("(",round(ci[i][1], sigdigits = 4), "-", round(ci[i][2], sigdigits = 4), ")"), halign = :left, family = "Palatino", pointsize = ps)
                annotate!(tp, 2.2, sty - i, t)
            end
        end

        l = @layout [a b{0.50w}]
        return plot(tp, p; layout = l, legend = false)
    else
        return plot!(p; kwargs...)
    end
end




forestplot([[0.2,1.2], [0.8,1.4], [0.6, 0.8], [1.2, 1.6]], 
       metric = [1.0, 1.2, 0.7, 1.3], source = ["A", "B", "C", "D"],
       summary = Dict(:ci =>[0.8, 1.1], :est => 0.95, :markershape => :rtriangle), logscale = false, refline=0.25, sourcelabel="coverage", metriclabel="estimate")


function afronden(x::Float64) 
    if isinteger(x)
        return string(Int(x))
    else
        y = round(x, sigdigits=3)
        if isinteger(y)
            return string(Int(y))
        else
            return y 
        end
    end 
end 
afronden(x::Int) = string(x)
afronden(x) = x 

results_afgerond_sorted = afronden.(results_sorted)

latexify(results_afgerond_sorted[results_afgerond_sorted.truenumind .== "2", [:truenumind, :coverage, :τC_low, :τC_median, :τC_high]], env=:table)

latexify(results_afgerond_sorted[results_afgerond_sorted.truenumind .== "1", [:truenumind, :coverage, :τA_low, :τA_median, :τA_high, :ϵ_low, :ϵ_median, :ϵ_high ]], env=:table)

latexify(results_afgerond_sorted[:, [:truenumind, :coverage, :τC_low, :τC_median, :τC_high, :τA_low, :τA_median, :τA_high]], env=:table)

latexify(results_afgerond_sorted[:, [:truenumind, :coverage, :ϵ_low, :ϵ_median, :ϵ_high]], env=:table)


@time @sync @distributed for (numberofindividualsinthesediment, coverage) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", coverage, "_MCMC.csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        println(zippedcsvfile, " does not exist.")
        continue
    end 
    plotfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "images", "subsampling", string("subsampling_number_of_individuals=", numberofindividualsinthesediment, "_coverage=", coverage, "_MCMC.png"))
    if isfile(plotfile) # plot is already made. 
        println("$plotfile does already exist. Skipping plotting.")
        continue 
    end 
    println("Working on ", plotfile)
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:truenumind, :coverage, :trueϵ, :chainid])
    chainarray = cat((Matrix(gdf[k][50_000:100_000, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
    chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
    p=plot(chain, show=false);
    savefig(p, plotfile);
end 

tableτA = Array{Vector{Float64}}(undef, 19, 8)
tableτC = Array{Vector{Float64}}(undef, 19, 8)
# table = zeros(19,8)

for (numind, coverage) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("knowntau=0.25_num_ind_", numind, "_cov_", coverage, "_standard_MCMC.csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:truenumind, :coverage, :trueϵ, :chainid])
    df = vcat((gdf[k][end÷2+1:10:end, :] for k in eachindex(gdf))...)
    tableτC[coverage-1, numind] = round.(quantile(df.:τCsample, [0.025, 0.5, 0.975]), digits=2)
    tableτA[coverage-1, numind] = round.(quantile(df.:τAsample, [0.025, 0.5, 0.975]), digits=2)
    #table[coverage-1, numind] = mean(df.τAsample)
end 

table = round.(table, digits=3)
latexify(table)

tabletuplesτC = similar(tableτC, Tuple{Float64, Float64, Float64})

for index in eachindex(tableτC, tabletuples)
    tabletuplesτC[index] = tuple(tableτC[index]...)
end 

tableτCstrings = string.(tableτC)
latexify(tableτCstrings)

tabletuplesτA = similar(tableτA, Tuple{Float64, Float64, Float64})

for index in eachindex(tableτA, tabletuplesτA)
    tabletuplesτA[index] = tuple(tableτA[index]...)
end 

tableτAstrings = string.(tabletuplesτA)
s = latexify(tableτAstrings[:, 5:8])

sp = split(s, '\n')
coverage = 2 
news = ""
for index in eachindex(sp)[4:end-4]
    news = news * string(coverage) * " & " * sp[index] * "\n\\hline\n"
    global coverage += 1 
end 
print(news)




tablecontained = similar(table, Bool)
for index in eachindex(table, tablecontained)
    tablecontained[index] = table[index][1] ≤ 0.25 ≤ table[index][3]
end



function printresults(anchorpop, ancientpop, truenumind, chainid, rate)
    output = string("anchor population ", anchorpop, " ancient population ", ancientpop, " true number of individuals ", truenumind, " chain ", chainid)
    output = output * string('\n', rate, "\n")
    println(output)
end 
# calculate acceptancerate
for (anchorpop, ancientpop, truenumind) in allcombinations
    csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("simprop_Madsmodel_ancient_pop_", ancientpop, "_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
    zippedcsvfile = csvfile * ".gz"
    if !isfile(zippedcsvfile)
        continue
    end 
    df = CSV.read(zippedcsvfile, DataFrame)
    gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :chainid])
    for df in gdf 
        rate = mean(df.:accepted[50_000:100_000])
        printresults(anchorpop, ancientpop, truenumind, df.:chainid[1], rate)
    end 
end 

anchorpop = "ACB"
truenumind = 2

csvfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "mcmcsamples", string("EmpericalData_anchor_pop_", anchorpop, "_true_num_ind_", truenumind, ".csv"))
zippedcsvfile = csvfile * ".gz"
    
zippedcsvfile = "data/mcmcsamples/storeliksEmpericalData_anchor_pop_MSL_true_num_ind_8.csv"
df = CSV.read(zippedcsvfile, DataFrame)

df = CSV.read(zippedcsvfile, DataFrame)#, header = [:anchorpop, :ancientpop, :trueϵ, :truenumind, :coverage, :nsample, :τCsample, :τAsample, :ϵsample, :accepted, :logjointprob, :chainid])

gdf = groupby(df, [:anchorpop, :ancientpop, :trueϵ, :truenumind, :chainid])

chainarray = cat((Matrix(gdf[k][5000:6000, [:nsample, :τCsample, :τAsample, :ϵsample ]]) for  k in 1:nchains)..., dims=3)
chain = Chains(chainarray, [:nsample, :τCsample, :τAsample, :ϵsample])
p=plot(chain)
