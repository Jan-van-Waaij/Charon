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



using GZip, # to work with zipped files,
      DataFrames, # to work with DataFrames,
      CSV, # for CSV files.
      StatsBase, # for the mean function
      CairoMakie # for plots 

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
end 

filename = "results.csv" 
CSV.write(filename, results)
cmd = `gzip -f $filename`
run(cmd)

gresults = groupby(results, :truenumind)

axvec = [(i, j) for i in 1:4 for j in 1:2]

parameters = [:n, :τC, :τA, :ϵ]
labels = [L"n", L"\tau_C", L"\tau_A", L"\epsilon"]
filenameparts = ["n", "tauC", "tauA", "eps"]
truevalue = [0.0, 0.25, 0.25, 0.003]

for parameterindex in eachindex(parameters, labels)
    f = Figure(; size = (1600, 1600))
    parameterlow = Symbol(parameters[parameterindex], :_low)
    parametermedium = Symbol(parameters[parameterindex], :_median)
    parameterhigh = Symbol(parameters[parameterindex], :_high) 
    @show parameterlow
    @show parametermedium
    @show parameterhigh
    for truenumind in 1:8
        @assert unique(gresults[truenumind].:truenumind) == [truenumind]
        ax = Axis(f[axvec[truenumind]...] , xlabel = "coverage", ylabel = labels[parameterindex], title = "Number of individuals: $truenumind")
        xs = gresults[truenumind].:coverage
        ys = gresults[truenumind][:, parametermedium]
        lowerrors = ys - gresults[truenumind][:, parameterlow]
        higherrors = gresults[truenumind][:, parameterhigh] - ys 
        errorbars!(ax, xs, ys, lowerrors, higherrors, whiskerwidth = 10, color = :red)
        scatter!(ax, xs, ys, markersize = 5, color = :black)
        if parameterindex == 1 
            hlines!(ax, truenumind, color = :blue )
        else 
            hlines!(ax, truevalue[parameterindex], color=:blue)
        end
    end 
    save(string("credibleintervals", filenameparts[parameterindex], ".png"), f)
end

for parameterindex in eachindex(parameters, labels, truevalue)
    f = Figure(; size = (1600, 1600))
    parameterlow = Symbol(parameters[parameterindex], :_low)
    parametermedium = Symbol(parameters[parameterindex], :_median)
    parameterhigh = Symbol(parameters[parameterindex], :_high) 
    @show parameterlow
    @show parametermedium
    @show parameterhigh
    for truenumind in 1:8
        @assert unique(gresults[truenumind].:truenumind) == [truenumind]
        ax = Axis(f[axvec[truenumind]...] , xlabel = "coverage", ylabel = labels[parameterindex], title = "Number of individuals: $truenumind")
        ys = gresults[truenumind].:coverage
        xs = gresults[truenumind][:, parametermedium]
        lowerrors = xs - gresults[truenumind][:, parameterlow]
        higherrors = gresults[truenumind][:, parameterhigh] - xs 
        errorbars!(ax, xs, ys, lowerrors, higherrors; whiskerwidth = 10, color = :red, direction = :x)
        scatter!(ax, xs, ys; markersize = 5, color = :black)
        if parameterindex == 1 
            vlines!(ax, truenumind, color = :blue )
        else 
            vlines!(ax, truevalue[parameterindex], color=:blue)
        end 
    end 
    save(string("credibleintervals", filenameparts[parameterindex], ".png"), f)
end


f = Figure(; size = (1600, 1600))
truenumind = 1 
ax = Axis(f[1,1] , xlabel = "coverage")
ys = gresults[truenumind].:coverage
xs = gresults[truenumind][:, :τC_median]
lowerrors = xs - gresults[truenumind][:, :τC_low]
higherrors = gresults[truenumind][:, :τC_high] - xs 
errorbars!(ax, xs, ys, lowerrors, higherrors; whiskerwidth = 10, color = :red, direction = :x)
scatter!(ax, xs, ys; markersize = 5, color = :black)
vlines!(ax, 0.25; ymin = minimum(ys), ymax = maximum(ys))


#=
maxcov = 10
fig = Figure()
ax = Axis(fig[1,1], yticks = 0.0:0.1:1.0, title = L"Credible sets for $\tau_A$")
#leg = Legend(fig[1,2])
#colbar = Colorbar(fig[1,2], limits=(1,8), colormap = cgrad(:viridis, 4, categorical = true))
bands = []
numindvec = [1, 4, 8]
for numind in numindvec
    topline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :τA_high]
    bottomline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :τA_low]
    coverages = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :coverage]
    
    push!(bands, band!(ax, coverages, bottomline, topline))
end 
Legend(fig[1,2], bands, string.(numindvec) .* [" individual"; fill(" individuals", 2)], title = "Number of individuals")
hlines!(ax, 0.25, linestyle = :dash )
fig 

save("data/images/crediblesets/crediblesetstauA.png", fig)


maxcov = 20
fig = Figure()
ax = Axis(fig[1,1],  title = L"Credible sets for $\tau_C$")
#leg = Legend(fig[1,2])
#colbar = Colorbar(fig[1,2], limits=(1,8), colormap = cgrad(:viridis, 4, categorical = true))
bands = []
numindvec = [1, 4, 8]
for numind in numindvec
    topline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :τC_high]
    bottomline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :τC_low]
    coverages = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :coverage]
    
    push!(bands, band!(ax, coverages, bottomline, topline))
end 
Legend(fig[1,2], bands, string.(numindvec) .* [" individual"; fill(" individuals", 2)], title = "Number of individuals")
hlines!(ax, 0.25, linestyle = :dash )
fig 

save("data/images/crediblesets/crediblesetstauC.png", fig)



maxcov = 10
fig = Figure(colormap = :viridis)
ax = Axis(fig[1,1],  title = L"Credible sets for $n$")
bands = []
numindvec = [1, 2, 3]
for numind in numindvec
    topline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :n_high]
    bottomline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :n_low]
    coverages = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :coverage]
    
    push!(bands, band!(ax, coverages, bottomline, topline))
end 
Legend(fig[1,2], bands, string.(numindvec) .* [" individual"; fill(" individuals", length(numindvec) - 1 )], title = "Number of individuals")
hlines!(ax, 1, linestyle = :dash )
hlines!(ax, 2, linestyle = :dash )
hlines!(ax, 3, linestyle = :dash )
fig 

save("data/images/crediblesets/crediblesetsn.png", fig)




maxcov = 10
fig = Figure()
ax = Axis(fig[1,1])
colbar = Colorbar(fig[1,2], limits=(1,8), colormap = cgrad(:viridis, 4, categorical = true))
for numind in [8, 4, 2, 1]
    topline = results[(results.:coverage .≤ maxcov), :τA_high]
    bottomline = results[(results.:coverage .≤ maxcov), :τA_low]
    coverages = results[ (results.:coverage .≤ maxcov), :coverage]
    numinds = results[(results.:coverage .≤ maxcov), :truenumind]

    @assert length(topline) == length(bottomline) == length(coverages)

    band!(ax, coverages, bottomline, topline, color = numinds)
end 
fig 



maxcov = 30
fig = Figure()
ax = Axis(fig[1,1],  title = L"Credible sets for $\epsilon$")
#leg = Legend(fig[1,2])
#colbar = Colorbar(fig[1,2], limits=(1,8), colormap = cgrad(:viridis, 4, categorical = true))
bands = []
numindvec = [1, 4, 8]
for numind in numindvec
    topline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :ϵ_high]
    bottomline = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :ϵ_low]
    coverages = results[(results.:truenumind .== numind) .& (results.:coverage .≤ maxcov), :coverage]
    
    push!(bands, band!(ax, coverages, bottomline, topline))
end 
Legend(fig[1,2], bands, string.(numindvec) .* [" individual"; fill(" individuals", 2)], title = "Number of individuals")
hlines!(ax, 0.003, linestyle = :dash )
fig 

save("data/images/crediblesets/crediblesetseps.png", fig)




fig = Figure()
ax = Axis(fig[1,1])
colors = [:red, :blue, :orange, :green]
for numind in 4:-1:1
    topline = df[df.numind .== numind, :high]
    bottomline = df[df.numind .== numind, :low]
    coverages = df[df.numind .== numind, :coverage]
    band!(ax, coverages, bottomline, topline, color = colors[numind])
end 
fig 



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

=# 