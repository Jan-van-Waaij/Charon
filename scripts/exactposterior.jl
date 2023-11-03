using PosteriorSampling, Data, CSV, DataFrames, Plots

frequencyfile = "/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqanchor/two_demes_11000_9000_70000_0_freqanchor.txt.gz"

R = 5 

readsfile = string("/home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads/two_demes_11000_9000_70000_0_4_", R, "_reads.txt.gz")

y = parsefreqfile(frequencyfile)

d = parsereadfile(readsfile)[1]

nrange = 1:8
τrange = 0.01:0.01:2.0
yrange = 0.001:0.001:0.999

# function allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, ϵrange::AbstractVector{Float64}, R::Int, d::Vector{Int}, y::Vector{Float64})
dic = allliks(nrange, τrange, τrange, R, d, y)

df = DataFrame(n = Int[], τC = Float64[], τA = Float64[], loglik = Float64[])

for key in keys(dic)
    push!(df, (key...,dic[key]))
end  

CSV.write("data/fullposterior/fullposteriorcoverage5firsttry.csv", df)


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
