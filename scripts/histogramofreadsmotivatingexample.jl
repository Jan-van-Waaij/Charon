using Distributions, Plots 

function simreads(n)
    distk = Binomial(2n, 0.3)
    k = rand(distk)
    distd = Binomial(100, k/(2n))
    return rand(distd)
end

ps = Vector{Any}(undef, 5)

for n in 1:5
    ds = [simreads(2^(n-1)) for _ in 1:10_000_000]
    ps[n] = histogram(ds, normalize= :pdf, bins=-0.5:100.5, xaxis=0:10:100, label = false, fill=true, title=n)
end 

p = plot(ps..., layout=(5))

nsamples = 10_000_000
nmax = 10 
ds = Matrix{Int}(undef, nsamples, nmax)
for n in 1:10
    for i in 1:nsamples
        ds[i, n] = simreads(n)
    end 
end 

histogram(ds, layout=(5,2), normalize= :pdf, bins=0.5:10.5, xaxis=1:10, label = false, fill=true,)


