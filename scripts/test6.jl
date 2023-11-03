using Likelihood, LinearAlgebra, Distributions, Dates, SpecialFunctions, StaticArrays
using Likelihood: makeQ, makeQꜜ, calcmatrix

import Base: exp

for n in 1:100
    Q = EigenExpansion(Matrix(makeQ(n)))
    Qꜜ = EigenExpansion(Matrix(makeQꜜ(n)))
    binomialcoefficients = binomial.(2n, 0:2n)
    for τC in 0.01:0.01:0.99
        for τA in 0.01:0.01:0.99
            M = calcmatrix(Q, Qꜜ, τC, τA)
            p = Vector{Float64}(undef, 2n+1) 
            for y in 0.01:0.01:0.99
                h = [y^k*(1-y)^(2n-k) for k in 0:2n]
                mul!(p, M, h)
                p .*= binomialcoefficients
                p .= max.(p, 0.0)
                @assert isprobvec(p)
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n ≤ 13

using Likelihood, LinearAlgebra, Distributions, StaticArrays
using Likelihood: calcmatrix, makeQ, makeQꜜ
function throughouttest()
    yrange = 0.00001:0.00001:0.99999
    τrange = 0.0001:0.0001:2.0
    for n in 1:13
        Q = EigenExpansion(makeQ(n))
        Qꜜ = EigenExpansion(makeQꜜ(n))
        binomialcoefficients = SVector{2n+1,Int}(binomial.(2n,0:2n))
        for τC in τrange
            for τA in τrange
                M = calcmatrix(Q, Qꜜ, τC, τA)
                p = MVector{2n+1,Float64}(undef) 
                for y in yrange 
                    h = SVector{2n+1,Float64}([y^k*(1-y)^(2n-k) for k in 0:2n])
                    mul!(p, M, h)
                    p .*= binomialcoefficients
                    p .= max.(p,0.0)
                    @assert isprobvec(p)
                end 
            end 
        end 
        println("finished with ", n)
    end # gaat goed tot en met n ≤ 13
end 

throughouttest()




function testEigenExpansion()
    for n in 1:13
        Q = EigenExpansion(Matrix(makeQ(n)))
        Qꜜ = EigenExpansion(Matrix(makeQꜜ(n)))
        binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
        for τC in 0.01:0.01:0.99
            for τA in 0.01:0.01:0.99
                M = calcmatrix(Q, Qꜜ, τC, τA)
                p = Vector{Float64}(undef, 2n+1) 
                for y in 0.01:0.01:0.99
                    h = [y^k*(1-y)^(2n-k) for k in 0:2n]
                    mul!(p, M, h)
                    p .*= binomialcoefficients
                end 
            end 
        end 
    end # gaat goed tot en met n ≤ 13
end # 52 sec

    
for n in 1:100
    Q = EigenExpansion(Matrix(makeQ(n)))
    Qꜜ = EigenExpansion(Matrix(makeQꜜ(n)))
    binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
    for τC in 0.01:0.01:0.99
        for τA in 0.01:0.01:0.99
            M = calcmatrix(Q, Qꜜ, τC, τA)
            p = Vector{Float64}(undef, 2n+1) 
            for y in 0.01:0.01:0.99
                h = [y^k*(1-y)^(2n-k) for k in 0:2n]
                mul!(p, M, h)
                p .*= binomialcoefficients
                @assert isprobvec(p)
            end 
        end 
    end 
end # gaat goed tot en met n ≤ 13

function testmyexp()
    for n in 1:13
        Q = eigen(Matrix(makeQ(n)))
        Qꜜ = eigen(Matrix(makeQꜜ(n)))
        binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
        for y in 0.01:0.01:0.99
            h = [y^k*(1-y)^(2n-k) for k in 0:2n]
            for τC in 0.01:0.01:0.99
                for τA in 0.01:0.01:0.99
                    p = myexp(Qꜜ, h, τC)
                    p = myexp(Q, p, τA)
                    p = binomialcoefficients .* p 
                end 
            end 
        end 
    end # gaat goed tot en met n ≤ 13
end # 453 sec

function testmyexp2()
    for n in 1:13
        Q = eigen(Matrix(makeQ(n)))
        Qꜜ = eigen(Matrix(makeQꜜ(n)))
        binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
        for y in 0.01:0.01:0.99
            h = [y^k*(1-y)^(2n-k) for k in 0:2n]
            for τC in 0.01:0.01:0.99
                for τA in 0.01:0.01:0.99
                    p = myexp2(Qꜜ, h, τC)
                    p = myexp2(Q, p, τA)
                    p = binomialcoefficients .* p 
                end 
            end 
        end 
    end # gaat goed tot en met n ≤ 13
end # 452 sec.

function testnaief()
    for n in 1:13
        Q = Matrix(makeQ(n))
        Qꜜ = Matrix(makeQꜜ(n))
        binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
        for y in 0.01:0.01:0.99
            h = [y^k*(1-y)^(2n-k) for k in 0:2n]
            for τC in 0.01:0.01:0.99
                for τA in 0.01:0.01:0.99
                    p = exp(τC*Qꜜ)*h
                    p = exp(τA*Q)*p
                    p = binomialcoefficients .* p  
                end 
            end 
        end 
    end # gaat goed tot en met n ≤ 13
end

@time testnaief() # 5859 sec. Super langzaam dus.

for n in 1:100
    Q = eigen(Matrix(makeQ(n)))
    Qꜜ = eigen(Matrix(makeQꜜ(n)))
    binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
    for y in 0.1:0.1:1.0
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = myexp2(Qꜜ, h, τC)
                p = myexp2(Q, p, τA)
                p = binomialcoefficients .* p 
                minp = minimum(p)
                sump = sum(p)
                @assert abs(sump-1.0) < 0.1 && minp > -0.1
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n\le 24

for n in 1:100
    Q = eigen(Matrix(makeQ(n)))
    Qꜜ = eigen(Matrix(makeQꜜ(n)))
    binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
    for y in 0.1:0.1:1.0
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = myexp(Qꜜ, h, τC)
                p = myexp(Q, p, τA)
                p = binomialcoefficients .* p 
                minp = minimum(p)
                sump = sum(p)
                @assert abs(sump-1.0) < 0.1 && minp > -0.1
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n\le 24

for n in 1:100
    Q = eigen(Matrix(makeQ(n)))
    Qꜜ = eigen(Matrix(makeQꜜ(n)))
    binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
    for y in 0.1:0.1:1.0
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = myexp(Qꜜ, h, τC)
                p = myexp(Q, p, τA)
                p = binomialcoefficients .* p 
                @assert isprobvec(p)
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n\le 13

for n in 1:100
    Q = eigen(Matrix(makeQ(n)))
    Qꜜ = eigen(Matrix(makeQꜜ(n)))
    binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
    for y in 0.1:0.1:1.0
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = myexp2(Qꜜ, h, τC)
                p = myexp2(Q, p, τA)
                p = binomialcoefficients .* p 
                @assert isprobvec(p)
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n\le 13

using Likelihood: makeQ, makeQꜜ
using LinearAlgebra, Distributions

for n in 1:100
    Q = Matrix(makeQ(n))
    Qꜜ = Matrix(makeQꜜ(n))
    binomialcoefficients = binomial.(2n,0:2n)
    for y in 0.1:0.1:1.0
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = exp(τC*Qꜜ)*h
                p = exp(τA*Q)*p
                p = binomialcoefficients .* p 
                minp = minimum(p)
                sump = sum(p)
                @assert isprobvec(p)
            end 
        end 
    end 
    println("finished with ", n)
end 


for n in 1:10
    Q = Matrix(makeQ(n))
    F = eigen(Q)
    println(maximum(abs.(Q - F.vectors*Diagonal(F.values)*inv(F.vectors)))) 
end 
            @assert normalexp(Qꜜ, x, τ) ≈ myexp(FQꜜ, x, τ)
function probderivedalleles(h::AbstractVector{Float64}, binomialcoefficients::AbstractVector{Float64}, Q::Eigen, Qꜜ::Eigen, τC::Float64, τA::Float64)
    a = exp(Diagonal(τA*Qꜜ.values))
end 


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
    n = (size(Q.P,1)-1)÷2
    for τC in τCrange
        for τA in τArange
            M = exp(Q, τA)*exp(Qꜜ, τC)
            for y in yrange 
                test(p, M, y)
                if abs(sum(p)-1.0)>0.01 || minimum(p) < -0.01
                    @show τC, τA, y, sum(p), minimum(p)
                end
            end 
        end 
    end 
end 

function overovertest(n, τrange, yrange)
    println("testing ", n)
    Q = NewEigenExpansion(makeQ(n))
    Qꜜ = NewEigenExpansion(makeQꜜ(n))
    p = MVector{2n+1,Float64}(undef)
    overtest(Q, Qꜜ, p, τrange, τrange, yrange)
end 

nrange = 19:100
yrange = 0.1:0.1:1.0 # vcat(0.001:0.001:0.009, 0.01:0.01:0.09, 0.1:0.1:0.9, 0.91:0.01:0.99, 0.991:0.001:0.999)
τrange = 0.1:0.1:1.0#vcat(0.001:0.001:0.003, 0.01:0.01:0.03, 0.1:0.1:0.9, 0.97:0.01:0.99, 0.999:0.001:1.002, 2.0, 10.0, 100.0)
for n in nrange
    overovertest(n, τrange, yrange)
end 
