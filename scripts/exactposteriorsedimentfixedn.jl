using PosteriorSampling, Data, CSV, DataFrames, GZip

function determinelocalmaxima(n, τCrange, τArange, ϵrange, dic)
    localmaxima = Tuple{Int, Float64, Float64, Float64}[]
    localmaximaindices = NTuple{3, Int64}[]
    localmaximavalues = Float64[]
    for τCindex in eachindex(τCrange)
        for τAindex in eachindex(τArange) 
            for ϵindex in eachindex(ϵrange)
                val = dic[(n, τCrange[τCindex], τArange[τAindex], ϵrange[ϵindex])]
                ismaximum = true
                for τCindex2 in max(1, τCindex-1):min(τCindex+1, length(τCrange))
                    for τAindex2 in max(1, τAindex-1):min(τAindex+1, length(τArange))
                        for ϵindex2 in max(1, ϵindex-1):min(ϵindex+1, length(ϵrange))
                            if (τCindex2, τAindex2, ϵindex2) != (τCindex, τAindex, ϵindex)
                                if dic[(n, τCrange[τCindex2], τArange[τAindex2], ϵrange[ϵindex2])] ≥ val 
                                    ismaximum = false 
                                end
                            end
                        end
                    end 
                end 
                if ismaximum
                    push!(localmaximaindices, (τCindex, τAindex, ϵindex))
                    push!(localmaxima, (n, τCrange[τCindex], τArange[τAindex], ϵrange[ϵindex]))
                    push!(localmaximavalues, val)
                end  
            end 
        end 
    end  
    return localmaxima, localmaximavalues
end 


vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles")
allfilesunfiltered = readdir(vcffilespath)
allfiles = String[]
for file in allfilesunfiltered 
    if file[end-6:end] != ".vcf.gz" # remove non VCF files.
        continue
    end 
    sp = split(file, '_')
    if sp[3] == "1000" # remove two_demes_1000_10_0.vcf.gz
        continue 
    end 
    push!(allfiles, file)
end 

fileswithoutextensionvec = [file[begin:end-7] for file in allfiles]



numberofindividualsinthesedimentvec = 1:8 #  [1, 2, 4, 8]  #  [1, 2, 4, 8, 16, 32] #, 16, 32, 64] 
nrange = 1:8 
τrange = 0.01:0.01:0.5
ϵrange = 0.0:0.0001:0.005

df = DataFrame(NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    trueτC = Float64[],
                    trueτA = Float64[],
                    trueϵ = Float64[],
                    truenumind = Int[],
                    n = Int[],
                    nummaxima = Int[],
                    localmaximavalues = Vector{Float64}[],
                    minlikelihood = Float64[],
                    maxlikelihood = Float64[],
                    localmaxima = Vector{Tuple{Int64, Float64, Float64, Float64}}[],
                    uniquemaximum = Bool[]
                    )

#df[:, [:NeAnchorPop, :NeAncientPop, :jointtimeyears, :trueτC, :trueτA, :trueϵ, :truenumind, :nummaxima]]

for filewithoutextension in fileswithoutextensionvec
    # filewithoutextension is of the form two_demes_Ne1_Ne2_jointimeyears_0
    sp = split(filewithoutextension, '_')
    NeAncientPop = parse(Int, sp[3])
    NeAnchorPop = parse(Int, sp[4])
    jointtimeyears = parse(Int, sp[5])
    generationtime = 25
    trueτA = jointtimeyears / (2*NeAncientPop*generationtime)
    trueτC = jointtimeyears / (2*NeAnchorPop*generationtime)
    frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", 
    "freqanchor", string(filewithoutextension, "_freqanchor.txt.gz"))
    gzopen(frequencyfile) do f 
        global freqs = parsefreqfile(f)
    end
    for truenumind in numberofindividualsinthesedimentvec
        readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", 
            string(filewithoutextension, '_', truenumind, '_', "1X_with_errors.txt.gz"))
        if !isfile(readsfile)
            println(readsfile, " does not exist.")
            continue
        end 
        gzopen(readsfile) do f 
            global derivedreadsvec, coveragevec = parsereadfile(f)
        end
        @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
        dic = allliks(nrange, τrange, τrange, ϵrange, coveragevec, derivedreadsvec, freqs)
        for n in nrange 
            localmaxima, localmaximavalues = determinelocalmaxima(n, τrange, τrange, ϵrange, dic)
            minlikelihood, maxlikelihood = extrema(values(dic))
            @show NeAnchorPop NeAncientPop jointtimeyears trueτC trueτA 0.003 truenumind n length(localmaxima) localmaximavalues minlikelihood maxlikelihood localmaxima 
            @assert length(localmaxima) == length(localmaximavalues)
            push!(df, (NeAnchorPop, NeAncientPop, jointtimeyears, trueτC, trueτA, 0.003, truenumind, n, length(localmaxima), localmaximavalues, minlikelihood, maxlikelihood, localmaxima, length(localmaxima) == 1))
        end 
    end 
end 

#=
df = DataFrame(NeAnchorPop = Int[],
                    NeAncientPop = Int[],
                    jointtimeyears = Int[],
                    trueτC = Float64[],
                    trueτA = Float64[],
                    trueϵ = Float64[],
                    truenumind = Int[],
                    nummaxima = Int[],
                    localmaximavalues = Vector{Float64}[],
                    localmaxima = Vector{Tuple{Int64, Float64, Float64, Float64}}[],
                    uniquemaximum = Bool[]
                    )
=#

#= saving for later use. 
CSV.write(data/fullposterior/fullposteriorsedimentdata.csv", df)
df2 = CSV.read(DataFrame, "data/fullposterior/fullposteriorsedimentdata.csv")
df2.localmaxima = eval.(Meta.parse.(df2.localmaxima))
=#



#=
#dic = allliks(nrange, τrange, τrange, R, d, y)
# (n, τC, τA, ϵ)
df = DataFrame(n = Int[], τC = Float64[], τA = Float64[], ϵ = Float64[], loglik = Float64[])

for key in keys(dic)
    push!(df, (key...,dic[key]))
end  

CSV.write("data/fullposterior/fullposteriorwitherrors.csv", df)

using CSV, DataFrames 

df = CSV.read("data/fullposterior/fullposteriorwitherrors.csv", DataFrame)

nrows = size(df,1)
dic = Dict((df[k,1], df[k,2], df[k,3], df[k,4]) => df[k,5] for k in 1:nrows)




df = CSV.read("data/fullposterior/fullposteriorcoverage5firsttry.csv", DataFrame)

nrows = size(df,1)
dic = Dict((df[k,1], df[k,2], df[k,3]) => df[k,4] for k in 1:nrows)

values = Array{Float64}(undef, length(τrange), length(τrange), length(nrange))
for i in eachindex(nrange)
    for j in eachindex(τrange)
        for k in eachindex(τrange)
            values[k, j, i] = dic[nrange[i], τrange[j], τrange[k]]
        end 
    end 
end

function islocalmaximum(x, i, j)
    xmax, ymax = size(x)
    if 1 < i < xmax
        if 1 < j < ymax 
            return x[i,j] > max(x[i-1,j-1], x[i-1,j], x[i-1,j+1], x[i,j-1], x[i,j+1], x[i+1,j-1], x[i+1,j], x[i+1,j+1])
        elseif j == 1
            return x[i,j] > max(x[i-1,j], x[i-1,j+1], x[i,j+1], x[i+1,j], x[i+1,j+1])
        elseif j == ymax 
            return x[i,j] > max(x[i-1,j-1], x[i-1,j], x[i,j-1], x[i+1,j-1], x[i+1,j])
        end 
    elseif i == 1 
        if 1 < j < ymax 
            return x[i,j] > max(x[i,j-1], x[i,j+1], x[i+1,j-1], x[i+1,j], x[i+1,j+1])
        elseif j == 1
            return x[i,j] > max(x[i,j+1], x[i+1,j], x[i+1,j+1])
        elseif j == ymax 
            return x[i,j] > max(x[i,j-1], x[i+1,j-1], x[i+1,j])
        end 
    elseif i == xmax 
        if 1 < j < ymax 
            return x[i,j] > max(x[i-1,j-1], x[i-1,j], x[i-1,j+1], x[i,j-1], x[i,j+1])
        elseif j == 1 
            return x[i,j] > max(x[i-1,j], x[i-1,j+1], x[i,j+1])
        elseif j == ymax 
            return x[i,j] > max(x[i-1,j-1], x[i-1,j], x[i,j-1])
        end 
    end 
    error("You forgot an option.")
end 

function countlocalmaxima(x)
    numlocmax = 0
    for i in 1:nrows
        for j in 1:ncols
            numlocmax += islocalmaximum(x, i, j)
        end 
    end 
    return numlocmax
end 

function islocalmaximum(x, i0, j0, k0)
    imax, jmax, kmax = size(x)
    islarger = trues(3,3,3)
    for i in max(i0-1,1):min(i0+1,imax)
        for j in max(j0-1,1):min(j0+1,jmax)
            for k in max(k0-1,1):min(k0+1,kmax)
                if (i, j, k) != (i0, j0, k0)
                    islarger[i-i0+2, j-j0+2, k-k0+2] = x[i0, j0, k0] > x[i, j, k]
                end 
            end 
        end 
    end
    return all(islarger) 
end  
    
function  countlocalmaxima(x)
    imax, jmax, kmax = size(x)
    numlocalmaxima = 0
    positions = Tuple{Int, Int, Int}[]
    for i in 1:imax
        for j in 1:jmax
            for k in 1:kmax
                if islocalmaximum(x, i, j, k)
                    numlocalmaxima += 1
                    push!(positions, (i,j,k))
                end  
            end 
        end 
    end 
    return numlocalmaxima, positions
end

nmax, maxsind = countlocalmaxima(values)
maxs = Vector{Tuple{Float64,Float64,Int}}(undef, nmax)
for index in 1:nmax   
    maxs[index] = τrange[maxsind[index][1]], τrange[maxsind[index][2]], nrange[maxsind[index][3]]
end 

#=
numlevels = 80
levels = range(extrema(values)..., length=numlevels)
p = similar(nrange, Any)
for n in eachindex(p)
    p[n] = contour(τrange, τrange, values[:, :, n], fill=true, levels=levels, title = string("n = ", n), cbar=false, aspect_ratio = :equal, xvals=τrange, yvals=τrange)
end
plot(p[1:8]..., layout=(4,2), xlabel = "τC", ylabel = "τA", size=(800,400))
=#

clims = extrema(values)
p = similar(nrange, Any)
for n in 1:7
    p[n] = heatmap(τrange, τrange, values[:, :, n], title = string("n = ", n), cbar=false, aspect_ratio = :equal, xlims=(0.0,2.0), ylims=(0.0,2.0), clims=clims)
end
p[8] = heatmap(τrange, τrange, values[:, :, 8], title = string("n = ", 8), aspect_ratio = :equal, xlims=(0.0,2.0), ylims=(0.0,2.0), clims=clims)
plot(p..., layout=(3,3), size=(900,900))
plot(p[1:8]..., layout=(4,2), size=(800,200)) #, xlabel = "τC", ylabel = "τA", size=(800,400))

gdf = groupby(df, :n)
=#