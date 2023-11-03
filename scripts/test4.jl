using Likelihood, LinearAlgebra, Distributions, Dates, SpecialFunctions
using Likelihood: makeQ, makeQꜜ

function custombinomial(n, k)
    try 
        return binomial(n,k)
    catch 
        return exp(logabsbinomial(n,k)[1])
    end 
end 

function test(p, M, y)   
    n = (size(M,1)-1)÷2     
    h = [y^k*(1-y)^(2n-k) for k in 0:2n]
    binomialcoefficients = custombinomial.(2n, 0:2n)
    mul!(p, M, h)
    p .*= binomialcoefficients
    return p 
end

function overtest(Q, Qꜜ, p, τCrange, τArange, yrange)
    n = (size(Q,1)-1)÷2
    for τC in τCrange
        for τA in τArange
            M = exp(τA*Q)*exp(τC*Qꜜ)
            for y in yrange 
                test(p, M, y)
                if !all(p .≥ 0.0) || !isprobvec(p) 
                    @show τC, τA, y, sum(p), minimum(p)
                end
            end 
        end 
    end 
end 

function overovertest(n, τrange, yrange)
    println("testing ", n)
    Q = Matrix(makeQ(n))
    Qꜜ = Matrix(makeQꜜ(n))
    p = Vector{Float64}(undef,2n+1)
    overtest(Q, Qꜜ, p, τrange, τrange, yrange)
end 

nrange = 34:100
yrange = vcat(0.001:0.001:0.009, 0.01:0.01:0.09, 0.1:0.1:0.9, 0.91:0.01:0.99, 0.991:0.001:0.999)
τrange = 0.1:0.1:100.0
for n in nrange
    overovertest(n, τrange, yrange)
end 

before = now()
Qs = [Matrix(makeQ(n)) for n in 1:100]
Qꜜs = [Matrix(makeQꜜ(n)) for n in 1:100] 
for turn in Base.OneTo(100_000)
    n = rand(1:100)
    τC = rand(0.1:0.1:100.0)
    τA = rand(0.1:0.1:100.0)
    M = exp(τA*Qs[n])*exp(τC*Qꜜs[n])
    @show turn
end # This took 4 hours and 23 minutes. 
after = now()
println(after - before)

println("finished")


for turn in Base.OneTo(100_000)
    n = rand(1:33)
    τC = rand(0.1:0.1:100.0)
    τA = rand(0.1:0.1:100.0)
    M = exp(τA*Qs[n])*exp(τC*Qꜜs[n])
    @show turn
end # This took 4 hours and 23 minutes. 
after = now()
println(after - before)

println("finished")

for n in 1:100
    Q = Matrix(makeQ(n))
    Qꜜ = Matrix(makeQꜜ(n))
    binomialcoefficients = exp.(logabsbinomial.(2n, k)[1] for k in 0:2n)
    for y in 0.1:0.1:0.9
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = exp(τC*Qꜜ)*h
                p = exp(τA*Q)*p 
                p = binomialcoefficients .* p
                minp = minimum(p)
                sump = sum(p)
                @assert abs(sump - 1.0) < 0.1 || minp > -0.1
            end 
        end 
    end 
    println("finished with ", n)
end 

# testing isprobvec gaat goed tot en met n ≤ 70. 

for n in 71:100
    Q = Matrix(makeQ(n))
    Qꜜ = Matrix(makeQꜜ(n))
    binomialcoefficients = exp.(logabsbinomial.(2n, k)[1] for k in 0:2n)
    for y in 0.1:0.1:0.9
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = exp(τC*Qꜜ)*h
                p = exp(τA*Q)*p 
                p = binomialcoefficients .* p
                @assert minimum(p) ≥ 0.0
                #minp = minimum(p)
                #sump = sum(p)
                #@assert abs(sump - 1.0) < 0.1 || minp > -0.1
            end 
        end 
    end 
    println("finished with ", n)
end 
