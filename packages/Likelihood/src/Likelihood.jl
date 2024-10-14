"""
A module to calculate the likelihood from read data.
"""
module Likelihood

# To do:
# test, write tests.
# adapt to unequal biomasses.

import Base: ==, exp, isapprox # redefine exp, == and isapprox (≈) for EigenExpansion types. This requires an import statement.

# external packages that are used in this module. 
using Distributions, # for probability distributions
      LinearAlgebra, # for Linear Algebra functions and data types.  
      StaticArrays # for small matrices this can be faster, as it stores the size of the matrix in the type. Can apply some 
      # optimisations for some matrices. 

export logprobderivedreads!, # the log probability of the read vector. A.k.a. the log-likelihood.
       EigenExpansion, # Data type to store a decomposition A = XDX^{-1}, which allows for fast matrix exponentiation.
       preparedata # prepare the data needed for logprobderivedreads!, for fixed n. 

# Test written
"""
    EigenExpansion{S, T}

Representation of a matrix M=P*D*P^{-1}, where D is a diagonal matrix and P an invertible matrix. P and P^{-1} have type S and the vector of diagonal elements of D has type T. The fields are P, D and Pinv. 
"""
struct EigenExpansion{S, T}
    P::S
    D::Diagonal{Float64, T}
    Pinv::S
    # inner constructor.
    function EigenExpansion(M::AbstractMatrix{Float64})
        F = eigen(M)
        P = F.vectors
        D = Diagonal(F.values)
        Pinv = inv(P)
        return new{typeof(P), typeof(F.values)}(P, D, Pinv)
    end
end 

# outer constructor
# Test written
"""
    EigenExpansion(M::Tridiagonal{Float64, Vector{Float64}})
    EigenExpansion(M::Matrix{Float64})

Try to find a decomposition M = P*D*P^(-1), where D is a diagonal matrix
and P is an invertible matrix. Note that exp(M) = P*exp(D)*P^(-1).
Fields:
F.P::Matrix{Float64}
F.D::Diagonal{Float64, Vector{Float64}}
F.Pinv::Matrix{Float64} 
"""
# This function can be problematic, because I don't always want an SMatrix. 
EigenExpansion(M::Tridiagonal{Float64, Vector{Float64}}) = EigenExpansion(SizedMatrix{size(M)...,Float64}(M))

# Test written. 
"""
    ==(E1::EigenExpansion, E2::EigenExpansion)

Test whether E1 and E2 have the same decomposition P*D*P^(-1), so whether E1.P == E2.P, E1.D == E2.D and E1.Pinv == E2.Pinv.
"""
==(E1::EigenExpansion, E2::EigenExpansion) = E1.P == E2.P && E1.D == E2.D && E1.Pinv == E2.Pinv

# test written 
"""
    isapprox(E1::EigenExpansion, E2::EigenExpansion; args...)
    E1 ≈ E2
    ≈(E1, E2; args...)

Test whether E1 and E2 have approximately the same decomposition P*D*P^(-1), so whether E1.P ≈ E2.P, E1.D ≈ E2.D and E1.Pinv ≈ E2.Pinv.
"""
isapprox(E1::EigenExpansion, E2::EigenExpansion, args...) = isapprox(E1.P, E2.P, args...) && isapprox(E1.D, E2.D, args...) && isapprox(E1.Pinv, E2.Pinv, args...)

# test written
"""
    exp(F::EigenExpansion, s::Float64)

Calculate e^(s*M) in an efficient way. 

If F=EigenExpansion(M), then exp(s*M)=exp(F,s)=F.P*exp(s*F.D)*F.Pinv, which makes 
use of the fast implementation of exp for diagonal matrices (it is just the exponetiation of the diagonal elements).     
"""
exp(F::EigenExpansion, s::Float64) =  F.P*exp(s*F.D)*F.Pinv
# Write a test function that tests whether typeof F.P is equal to the type of the end result.  DONE

# Test written.
"""
    makeQ(n::Int)

Calculate the (2n+1)x(2n+1) Q tridiagonal matrix as in Schraiber 2018. 
"""
function makeQ(n::Int)
    lowerdiagonal = 0.5*(1.0:2n).*(0.0:2n-1.0)
    diagonal = -(0.0:2n).*(2n:-1.0:0.0)
    upperdiagonal = 0.5*(2n:-1.0:1.0).*(2n-1:-1.0:0.0)
    return Tridiagonal(lowerdiagonal, diagonal, upperdiagonal)
end 

# Test written.
"""
    makeQꜜ(n::Int)

Calculate the (2n+1)x(2n+1) Qꜜ tridiagonal matrix as in Schraiber 2018.  
"""
function makeQꜜ(n::Int)
    lowerdiagonal = 0.5*(1.0:2n).*(0.0:2n-1.0)
    diagonal = -(0.0:2n).*(2n+1.0:-1.0:1.0)
    upperdiagonal = 0.5*(2n+1:-1.0:2.0).*(2n:-1.0:1.0)
    return Tridiagonal(lowerdiagonal, diagonal, upperdiagonal)
end

# Tests written. 
"""
    calcmatrix(Q::EigenExpansion, Qꜜ::EigenExpansion, τC::Float64, τA::Float64)

Calculate the matrix `e^{Q*τ_A}e^{Qꜜτ_C}` as in Schraiber. Note that in general, for matrices, exp(A+B) ≠ exp(A)exp(B) (equality holds if and only if A and B commute, which is not the case for Q and Qꜜ). 
"""
calcmatrix(Q::EigenExpansion, Qꜜ::EigenExpansion, τC::Float64, τA::Float64) =  exp(Q, τA)*exp(Qꜜ, τC)

# Test written
#=
"""
    probderivedalleles!(p::AbstractVector{Float64}, h::AbstractVector{Float64}, M::AbstractMatrix{Float64}, binomialcoefficients::AbstractVector{Int})

Updates the vector p in place, the other arguments remain unchanged. The vector p is of size 2n+1. Entry p[k+1] is the probability of k alleles, k=0,...,2n (remember Julia indexing starts at 1). `M=e^{Q*τ_A}e^{Qꜜτ_C}` is the matrix calculated with calcmatrix.
h is a (2n+1) x length(y)-matrix, where h[k+1, j] = y[j]^k*(1-y[j])^{2n-k}, k = 0,...,2n. binomialcoefficients is a vector with the k+1-element equal to the binomial coefficient 2n over k, k=0,...,2n.

The programmer is responsible for providing the correct data!
"""
function probderivedalleles!(p::AbstractVector{Float64}, h::AbstractVector{Float64}, M::AbstractMatrix{Float64}, binomialcoefficients::AbstractVector{Int})
    # mul! is 60% faster. 
    mul!(p, M, h) # equivalent to p = M*h, but in place, without allocating new vectors.
    p .*= binomialcoefficients # in place computation of p = p*binomialcoefficients
    p .= max.(0.0, p)
    @assert isprobvec(p)
    return p  
end
=#
# Het omslagpunt ligt bij length(x) ≤ 23 (dus n=11), voor length(x)≤23 is het sneller met StaticArrays,
# voor length(x)≥24 is een gewone array sneller (dus voor n≥12). 
#
# Maar wat nu als M relatief klein is, en p, h en binomialcoefficients zijn views uit mega-vectors?
#  Met een vector van grootte 1000 ligt het omslagpunt bij 20 ongeveer. Het resultaat is erg beetje wisselend.
# Het lijkt erop dat SizedArrays is een goede compromis. Als de length klein is, dan doet die het ietsjes 
# slechter dan StaticArrays, maar niet heel veel slechter. Daarentegen als n groot is, dan doet die het 
# maar ietsjes slechter dan normal arrays, maar aanzienlijk beter dan static arrays. 

# HIER VERDER MET TESTEN SCHRIJVEN!

#=
"""
    makebinomialprobs(n::Int, R::Int)

Returns an (2n+1)x(R+1) matrix of type SMatrix M where entry M[k+1, d+1] contains the probability of d derived reads out of R, given k derived reads out of 2n, which is P(X=d), with X binomially distributed with parameters R and k/2n.  
"""
function makebinomialprobs(n::Int, R::Int)
    # binomial distributions with R experiments and success probabilities 0, 1/(2n),..., 1.
    binomialdistributions = Binomial.(R, range(0.0, 1.0, length = 2n+1)) # with Distributions package. 
    binomialprobs = Matrix{Float64}(undef, 2n+1, R+1)
    for d in 0:R 
        for k in 1:2n+1
            binomialprobs[k, d+1] = pdf(binomialdistributions[k], d) 
        end 
    end 
    return SMatrix{2n+1, R+1, Float64}(binomialprobs) # type unstable, because the output type depends on the input n and R. But does not matter, we only use it once.
end 

"""
    makebinomialprobs(n::Int, R::Vector{Int})

Make a vector of length maximum(R) where the r-th element contains the matrix generated by makebinomialprobs(n, r).
"""
function makebinomialprobs(n::Int, R::Vector{Int})
    Rmax = maximum(R)
    binomialprobs = [makebinomialprobs(n, r) for r in 1:Rmax]
    return SVector{Rmax, SArray{S, Float64, 2, L} where {S<:Tuple, L}}(binomialprobs)
end 
=#

"""
    preparedata(n::Int, y::AbstractVector{<:Real})

Returns a named tuple (Q=Q, Qꜜ=Qꜜ, hvals=hvals, binomialcoefficients=binomialcoefficients), where 
Q and Qꜜ are EigenExpansions of Q and Qꜜ (as in Schraiber 2018), respectively. 
hvals is a (2n+1) x length(y) matrix of type SizedMatrix. SizedMatrix is used rather than SMatrix, because this matrix is so large, and 
SMatrix is slow for large matrices, here hvals[k+1,i] = y[i]^k*(1-y[i])^(2n-k).
binomialcoefficients is a vector (of type SVector) containing the binomial coefficients 2n over k, k=0,...,2n.  
"""
function preparedata(n::Int, y::AbstractVector{<:Real})
    # I can define an if statement and give Matrix(makeQ(n)) for large n.
    Q = EigenExpansion(makeQ(n))
    Qꜜ = EigenExpansion(makeQꜜ(n))
    binomialcoefficients = SVector{2n+1, Int}(binomial.(2n, 0:2n))
    hvals = SizedMatrix{2n+1, length(y), Float64}(undef)
    for j in eachindex(y)
        for k in 0:2n 
            hvals[k+1, j] = y[j]^k*(1.0-y[j])^(2n-k)
        end 
    end
    return (Q=Q, Qꜜ=Qꜜ, hvals=hvals, binomialcoefficients=binomialcoefficients)
end

## With errors and constant coverage.
#=
"""
    logprobderivedreads!(p::AbstractVector{Float64}, q::AbstractArray{Float64}, reads::AbstractVector{Int}, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, R::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})

Calculates the log-probability of the data (a.k.a. log-likelihood, data is d) given the parameters and the prior data. It is implicitly assumed that Q and Qꜜ are (2n+1)x(2n+1)-matrices (or their eigenvalue decomposition), that p is a 2n+1 vector, binomialcoefficients is a 2n+1-vector of binomial coefficients, with k+1 element equal to binomial coefficient 2n over k, k = 0,...,2n. binomialprobs is a (2n+1)x(R+1)-matrix where the (k+1, d+1) element is equal to the probability of d success out of R experiments with probability of success k/(2n), where k = 0,...,2n, and d = 0,...,R. hvals is a (2n+1)xlength(y) matrix, where y are all frequencies, where the (k+1,j)-th element is equal to y[j]^k*(1.0-y[j])^(2n-k). So you should initialise 
Q, Qꜜ, hvals, binomialcoefficients, p, binomialprobs = preparedatafixedn(n, R, y)
"""
function logprobderivedreads!(p::AbstractVector{Float64}, q::AbstractArray{Float64}, reads::AbstractVector{Int}, R::Int, Rmax::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
    if τC ≤ 0 || τA ≤ 0 || ϵ < 0 || ϵ > 1
        return -Inf 
    else
        M = calcmatrix(Q, Qꜜ, τC, τA)
        # update q. q[k+1, d+1] is the probability of d derived alleles given k derived alleles (out of 2n)
        # and error probability ϵ.
        for k in 0:2n
            dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
            for d in 0:R
                q[k+1, d+1] = pdf(dist, d)
            end 
        end 
        loglik = 0.0
        for index in eachindex(reads, counts)
            # update p
            probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
            loglik +=  counts[index]*log(dot(p, view(q, :, reads[index]+1))) 
        end
        return loglik
    end 
end # Constant coverage and with errors. 

## Implement likelihood for non-constant coverages and with errors.
function logprobderivedreads!(p::AbstractVector{Float64}, q::AbstractVector{<:AbstractMatrix{Float64}}, reads::AbstractVector{Int}, coverages::AbstractVector{Int}, maxcoverage::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
    if τC ≤ 0 || τA ≤ 0 || ϵ < 0 || ϵ > 1
        return -Inf
    else
        M = calcmatrix(Q, Qꜜ, τC, τA)
        for R in 1:maxcoverage
            for k in 0:2n
                dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
                for d in 0:R
                    q[R][k+1, d+1] = pdf(dist, d)
                end 
            end 
        end 
        loglik = 0.0
        for index in eachindex(reads, coverages, counts)
            # update p
            probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
            # p[k+1] now contains the probability of k derived alleles at locus index,
            # with n individuals and model parameters τC and τA and anchor frequency y[index].
            loglik +=  counts[index]*log(dot(p, view(q[coverages[index]], :, reads[index] + 1)))
        end
        return loglik
    end 
end
## End implementation voor verschillende coverages.

function logprobderivedreads!(p::AbstractVector{Float64}, q::Dict{Tuple{Int64, Int64}, Vector{Float64}}, reads::AbstractVector{Int}, coverages::AbstractVector{Int}, maxcoverage::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
    if τC ≤ 0 || τA ≤ 0 || ϵ < 0 || ϵ > 1
        return -Inf
    else
        M = calcmatrix(Q, Qꜜ, τC, τA) 
        for R in unique(coverages)
            for k in 0:2n
                dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
                for d in 0:R
                    if haskey(q, (d, R))
                        q[d, R][k+1] = pdf(dist, d)
                    end 
                end 
            end 
        end 
        loglik = 0.0
        for index in eachindex(reads, coverages, counts)
            # update p
            probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
            # p[k+1] now contains the probability of k derived alleles at locus index,
            # with n individuals and model parameters τC and τA and anchor frequency y[index].
            loglik += counts[index]*log(dot(p, q[reads[index], coverages[index]]))
        end
        return loglik
    end 
end
## End implementation voor verschillende coverages.
=#

#= 

=#
#=
function probderivedallelesreminder!(p::AbstractVector{Float64}, h::AbstractVector{Float64}, M::AbstractMatrix{Float64}, binomialcoefficients::AbstractVector{Int})
    # mul! is 60% faster. 
    mul!(p, M, h) # equivalent to p = M*h, but in place, without allocating new vectors.
    p .*= binomialcoefficients # in place computation of p = p*binomialcoefficients
    p .= max.(0.0, p)
    @assert isprobvec(p)
    return p  
end
=#


"""
    logprobderivedreads!(q::AbstractArray{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})

Calculate the log probability of the derived reads, which is the log-likelihood. This function will mutate q. q[k+1,d+1] contains the probability of seeing d derived reads out of R, given k derived alleles out of 2n, d = 0,...,R and k=0,...,2n.  
"""
function logprobderivedreads!(q::AbstractArray{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
    if τC < 0 || τA < 0 || ϵ < 0 || ϵ > 1
        return -Inf 
    else
        M = calcmatrix(Q, Qꜜ, τC, τA)
        # update q. q[k+1, d+1] is the probability of d derived alleles given k derived alleles (out of 2n)
        # and error probability ϵ.
        for k in 0:2n
            dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
            for d in 0:R
                q[k+1, d+1] = binomialcoefficients[k+1] * pdf(dist, d)
            end
        end 
        loglik = 0.0
        for index in eachindex(reads, counts)
            # update p
            loglik += counts[index] * log(max(0.0, dot(view(q, :, reads[index]+1), M, view(hvals, :, index))))
            #probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
            #loglik +=  counts[index]*log(dot(p, view(q, :, reads[index]+1))) 
        end
        return loglik
    end 
end # Constant coverage and with errors. 

"""
    logprobderivedreads!(q::Dict{Tuple{Int64, Int64}, Vector{Float64}}, reads::AbstractVector{Int}, coverages::AbstractVector{Int}, uniquecoverages::AbstractVector{Int}, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})

Same as before, only with non-constant coverages. Mutates q, which is dictionary. Key (d, R), if it exists, has an associate vector as value, where the k+1-th element is the probability of having d derived reads out of R, given k derived alleles out of 2n, d=0,...,R and k=0,...,2n. 
"""
function logprobderivedreads!(q::Dict{Tuple{Int64, Int64}, Vector{Float64}}, reads::AbstractVector{Int}, coverages::AbstractVector{Int}, uniquecoverages::AbstractVector{Int}, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
    if τC < 0 || τA < 0 || ϵ < 0 || ϵ > 1
        return -Inf
    else
        M = calcmatrix(Q, Qꜜ, τC, τA) 
        for R in uniquecoverages
            for k in 0:2n
                dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
                for d in 0:R
                    if haskey(q, (d, R))
                        q[d, R][k+1] = binomialcoefficients[k+1]*pdf(dist, d)
                    end 
                end 
            end 
        end 
        loglik = 0.0
        for index in eachindex(reads, coverages, counts)
            loglik += counts[index]*log(max(0.0, dot(q[reads[index], coverages[index]], M, view(hvals, :, index))))
        end
        return loglik
    end
end


###
### Unequal contributions
###


"""
    logprobderivedreads!(q::AbstractArray{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})

Calculate the log probability of the derived reads, which is the log-likelihood. This function will mutate q. q[k+1,d+1] contains the probability of seeing d derived reads out of R, given k derived alleles out of 2n, d = 0,...,R and k=0,...,2n.  
"""
function logprobderivedreadsunevencontributions!(probderreadgivenk::AbstractArray{Float64}, reads::AbstractVector{Int}, coverage::Int, uniqueCoverage::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, numind::Int, proportions::Vector{Float64}, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
    if τC < 0 || τA < 0 || ϵ < 0 || ϵ > 1
        return -Inf 
    else
        M = calcmatrix(Q, Qꜜ, τC, τA)
        # update q. q[k+1, d+1] is the probability of d derived alleles given k derived alleles (out of 2n)
        # and error probability ϵ.
        proportionslong = repeat(proportions, inner=2)/2
        for k in 0:2numind
            dist = Binomial(coverage, k*(1.0-ϵ)/(2numind) + (2numind-k)*ϵ/(2numind))
            for d in 0:coverage
                probderreadgivenk[k+1, d+1] = binomialcoefficients[k+1] * pdf(dist, d)
            end
        end 
        loglik = 0.0
        for index in eachindex(reads, counts)
            # update p
            loglik += counts[index] * log(max(0.0, dot(view(probderreadgivenk, :, reads[index]+1), M, view(hvals, :, index))))
            #probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
            #loglik +=  counts[index]*log(dot(p, view(q, :, reads[index]+1))) 
        end
        return loglik
    end 
end # Constant coverage and with errors. 

function generateall2nsequencessummingtok(n::Int, k::Int)
    allsequences = Vector{Bool}[]
    for m in 0:2^(2n)-1
        if count_ones(m) == k
            push!(allsequences, digits(Bool, m, base=2, pad=2n))
        end 
    end 
    return allsequences
end 


end # module Likelihood