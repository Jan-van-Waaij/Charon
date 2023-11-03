using Likelihood, LinearAlgebra, Distributions, Dates, SpecialFunctions, StaticArrays
using Likelihood: makeQ, makeQꜜ

import Base: exp

struct NewEigenExpansion{S, T}
    P::S
    D::Diagonal{Float64, T}
    Pinv::S
    # inner constructor.
    function NewEigenExpansion(M::AbstractMatrix{Float64})
        m = maximum(abs.(M))
        F = eigen(M/m)
        P = F.vectors
        D = Diagonal(m*F.values)
        Pinv = inv(P)
        return new{typeof(P), typeof(F.values)}(P, D, Pinv)
    end
end 

# This function can be problematic, because I don't always want an SMatrix. 
NewEigenExpansion(M::Tridiagonal{Float64, Vector{Float64}}) = NewEigenExpansion(Matrix(M))

# test written
"""
    exp(F::EigenExpansion, s::Float64)

If F=EigenExpansion(M), then exp(s*M)=exp(F,s)=F.P*exp(s*F.D)*F.Pinv, which makes 
use of the fast implementation of exp for diagonal matrices.     
"""
exp(F::NewEigenExpansion, s::Float64) =  F.P*exp(s*F.D)*F.Pinv
# Write a test function that tests whether typeof F.P is equal to the type of the end result.  DONE

function myexp(A::Eigen, x::AbstractVector{Float64}, τ::Float64)
    # Say, A = X^{-1}DX, where D is a diagonal matrix, and X and invertible matrix.
    # So X is the matrix with the eigenvectors as columns, and D is the diagonal matrix
    # where the d-th diagonal element is the eigenvalue corresponding to the eigenvector 
    # at the d-th column. 
    a = A.vectors \ x 
    b = A.vectors * exp(Diagonal(τ*A.values)) * a 
    return b 
end 

function myexp2(A::Eigen, x::AbstractVector{Float64}, τ::Float64)
    # Say, A = X^{-1}DX, where D is a diagonal matrix, and X and invertible matrix.
    # So X is the matrix with the eigenvectors as columns, and D is the diagonal matrix
    # where the d-th diagonal element is the eigenvalue corresponding to the eigenvector 
    # at the d-th column. 
    a = lu(A.vectors) \ x 
    b = A.vectors * exp(Diagonal(τ*A.values)) * a 
    return b 
end 

function myexp3!(A::Eigen, F::LU, x::AbstractVector{Float64}, τ::Float64)
    # Say, A = X^{-1}DX, where D is a diagonal matrix, and X and invertible matrix.
    # So X is the matrix with the eigenvectors as columns, and D is the diagonal matrix
    # where the d-th diagonal element is the eigenvalue corresponding to the eigenvector 
    # at the d-th column. 
    a = lu(A.vectors) \ x 
    b = A.vectors * exp(Diagonal(τ*A.values)) * a 
    return b 
end 


normalexp(A, x, τ) = exp(τ*A)*x

@time for n in 1:10
    Q = Matrix(makeQ(n))
    Qꜜ = Matrix(makeQꜜ(n))
    for y in 0.01:0.01:1.0
        x = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τ in 0.01:0.01:1.0
            normalexp(Q, x, τ)
            normalexp(Qꜜ, x, τ)
        end 
    end 
end 

@time for n in 1:10
    Q = Matrix(makeQ(n))
    FQ = eigen(Q)
    Qꜜ = Matrix(makeQꜜ(n))
    FQꜜ = eigen(Qꜜ)
    for y in 0.01:0.01:1.0
        x = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τ in 0.01:0.01:1.0
            myexp(FQ, x, τ)
            myexp(FQꜜ, x, τ)
        end 
    end 
end 

for n in 1:100
    Q = eigen(Matrix(makeQ(n)))
    Qꜜ = eigen(Matrix(makeQꜜ(n)))
    binomialcoefficients = exp.(logabsbinomial(2n,k)[1] for k in 0:2n)
    for y in 0.1:0.1:0.9
        h = [y^k*(1-y)^(2n-k) for k in 0:2n]
        for τC in 0.1:0.1:1.0
            for τA in 0.1:0.1:1.0
                p = myexp(Qꜜ, h, τC)
                p = myexp(Q, p, τA)
                p = binomialcoefficients .* p 
                minp = minimum(p)
                sump = sum(p)
                @assert abs(sump-1.0) < 0.1 || minp > -0.1
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n\le 25. 

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
                @assert abs(sump-1.0) < 0.1 || minp > -0.1
            end 
        end 
    end 
    println("finished with ", n)
end # gaat goed tot en met n\le 25

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
