"""
    module PackageToBeShipped

This package provides functions to sample from the posterior. 
"""
module PackageToBeShipped

import Base: ==, exp # redefine exp and == for EigenExpansion types. This requires an import statement.

# external packages that I use. 
using Distributions # for probability distributions
using StatsBase # for the countmap function
using LinearAlgebra # Needed for the identity matrix I. 
using StaticArrays # for SizedVector and MVector, they are generally faster.
using Base.Threads # for multithreading.
using DataFrames # for processing Dice files. 
using CSV # for processing CSV files. 
      
# functions that you can use when you type 'using PosteriorSampling'. 
export MCMCsampler # More efficient mcmc sampler for low coverage with simulated proposal, at the cost of a longer run time. 

## Likelihood functions, from the former Likelihood package 
# Tests written
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
# Tests written
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
EigenExpansion(M::Tridiagonal{Float64, Vector{Float64}}) = EigenExpansion(SizedMatrix{size(M)...,Float64}(M))

# Tests written. 
"""
    ==(E1::EigenExpansion, E2::EigenExpansion)

Test whether E1 and E2 have the same decomposition P*D*P^(-1), so whether E1.P == E2.P, E1.D == E2.D and E1.Pinv == E2.Pinv.
"""
==(E1::EigenExpansion, E2::EigenExpansion) = E1.P == E2.P && E1.D == E2.D && E1.Pinv == E2.Pinv

# tests written
"""
    exp(F::EigenExpansion, s::Float64)

Calculate e^(s*M) in an efficient way. 

If F=EigenExpansion(M), then exp(s*M)=exp(F,s)=F.P*exp(s*F.D)*F.Pinv, which makes 
use of the fast implementation of exp for diagonal matrices (it is just the exponetiation of the diagonal elements).     
"""
exp(F::EigenExpansion, s::Float64) =  F.P*exp(s*F.D)*F.Pinv

# Tests written.
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

# Tests written.
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

# Tests written 
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

# tests written 
"""
    logprobderivedreads!(q::AbstractArray{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})

Calculate the log probability of the derived reads, which is the log-likelihood. This function will mutate q. q[k+1,d+1] contains the probability of seeing d derived reads out of R, given k derived alleles out of 2n, d = 0,...,R and k=0,...,2n.  
"""
function logprobderivedreads!(q::AbstractMatrix{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
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
        end
        return loglik
    end 
end # Constant coverage and with errors. 


# tests written 
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


## End likelihood functions. 

function filtervectorsandapplycountmap(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}, allowedindices::AbstractVector{Bool})
    y, d, R = y[allowedindices], d[allowedindices], R[allowedindices]
    tuples = zip(d, R, y)
    cm = countmap(tuples)
    d = [key[1] for key in keys(cm)]
    R = [key[2] for key in keys(cm)]
    y = [key[3] for key in keys(cm)]
    counts = collect(values(cm))
    return y, d, R, counts
end 

function filtervectorsandapplycountmap(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::Integer, counts::AbstractVector{<:Integer}, allowedindices::AbstractVector{Bool})
    y, d, counts = y[allowedindices], d[allowedindices], counts[allowedindices]
    return y, d, R, counts
end 

function filtervectorsandapplycountmap(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, allowedindices::AbstractVector{Bool})
    y, d, R, counts = y[allowedindices], d[allowedindices], R[allowedindices], counts[allowedindices]
    return y, d, R, counts
end 

function makesizedvectors(d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    d = SizedVector{length(d), Int}(d)
    R = SizedVector{length(R), Int}(R)
    return d, R 
end 

function makesizedvectors(d::AbstractVector{<:Integer}, R::Integer)
    d = SizedVector{length(d), Int}(d)
    return d, R 
end 

function makeqfixedn(n::Integer, d::AbstractVector{<:Integer}, R::Integer)
    return MMatrix{2n+1, R+1, Float64}(undef)
end 

function makeqfixedn(n::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    uniquedandRs = unique(zip(d, R))
    return Dict((d, R) => Vector{Float64}(undef, 2n+1) for (d, R) in uniquedandRs)
end 

function makeq(nmax::Integer, d::AbstractVector{<:Integer}, R::Integer)
    return [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
end 

function makeq(nmax::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    uniquedandRs = unique(zip(d, R))
    return [Dict((d, R) => Vector{Float64}(undef, 2n+1) for (d, R) in uniquedandRs) for n in 1:nmax]
end 

checklength(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}) = length(y) == length(d) == length(R)

checklength(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::Integer) = length(y) == length(d)

checklength(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::Integer, counts::AbstractVector{<:Integer}) = length(y) == length(d) == length(counts)

checklength(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}) = length(y) == length(d) == length(counts) == length(R)

function updateq!(q::AbstractMatrix{<:Real}, binomialcoefficients::AbstractVector{<:Integer}, uniquecoverages::Integer, n::Integer, ϵ::Real)
    for k in 0:2n
        dist = Binomial(uniquecoverages, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
        for d in 0:uniquecoverages
            q[k+1, d+1] = binomialcoefficients[k+1] * pdf(dist, d)
        end 
    end              
end 

function updateq!(q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, binomialcoefficients::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, n::Integer, ϵ::Real)
    for R in uniquecoverages
        for k in 0:2n
            dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
            for d in 0:R
                if haskey(q, (d, R))
                    q[d, R][k+1] = binomialcoefficients[k+1] * pdf(dist, d)
                end 
            end
        end     
    end          
end 

calcloglik(counts::AbstractVector{<:Integer}, q::AbstractMatrix{<:Real}, reads::AbstractVector{<:Integer}, coverages::Union{Integer, AbstractVector{<:Integer}}, M::AbstractMatrix{<:Real}, hvals::AbstractMatrix{<:Real}, index::Integer) = counts[index] * log(max(0.0, dot(view(q, :, reads[index]+1), M, view(hvals, :, index))))

calcloglik(counts::AbstractVector{<:Integer}, q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, M::Matrix{<:Real}, hvals::Matrix{<:Real}, index::Integer) = counts[index]*log(max(0.0, dot(q[reads[index], coverages[index]], M, view(hvals, :, index))))

myunique(R::Integer) = R 

myunique(R::AbstractVector{<:Integer}) = unique(R)

function readcsvfile(filename::AbstractString)
    dftest = CSV.read(filename, DataFrame, header = false, limit=1)  
    if typeof(dftest[1,1]) <: AbstractString
        df = CSV.read(filename, DataFrame, header = true)
    else
        df = CSV.read(filename, DataFrame, header = false)
    end 
    return df 
end 

##
##
## fixed R, error model, fixed n. 
##
##
"""
    exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::Union{AbstractVector{<:Integer}, Integer}, uniquecoverages::Union{AbstractVector{<:Integer}, Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on τC, τA and ϵ, the true number of individuals is assumed to be known. Instead of using the recommended scaling factor 1.4161, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress.

If not specified, the initial values are choosen from the prior.

This function returns a vector of length nchains. Each element contains the results of one chain; it is a tuple of five vectors: a vector containing the a vector containing the τC-values, a vector containing the τA-values, a vector containing the ϵ values, a vector telling whether a certain proposal was accepted or not and a vector with unnormalised posterior log-densities. So for example 

```julia
results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, R, d, y, messages)
τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[1]
```

gives the results of the first chain.

Here we assume that the true number of individuals in the sediment is known and that the coverage is constant.
"""
function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::Union{AbstractVector{<:Integer}, Integer}, uniquecoverages::Union{AbstractVector{<:Integer}, Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer)
    # Check conditions.
    all(0 .≤ ϵrange .< 0.5) || throw(DomainError("ϵ should be in [0.0,0.5)."))
    length(nrange) > 0 || throw(ArgumentError("nrange should have at least one element."))
    length(τCrange) > 0 || throw(ArgumentError("τCrange should have at least one element."))
    length(τArange) > 0 || throw(ArgumentError("τArange should have at least one element."))
    length(ϵrange) > 0 || throw(ArgumentError("ϵrange should have at least one element."))
    
    nrows = length(nrange)*length(τCrange)*length(τArange)*length(ϵrange)
    ns = Vector{Int}(undef, nrows)
    τCs = Vector{Float64}(undef, nrows)
    τAs = Vector{Float64}(undef, nrows)
    ϵs = Vector{Float64}(undef, nrows)
    logliks = Vector{Float64}(undef, nrows)

    index = 1
    for n in nrange
        Q, Qꜜ, hvals, binomialcoefficients = preparedata(n, y)
        q = makeqfixedn(n, reads, coverages)
        M = Matrix{Float64}(undef, 2n+1, 2n+1)
        for τC in τCrange 
            expτCQꜜ = exp(Qꜜ, τC)
            for τA in τArange
                expτAQ = exp(Q, τA)    
                mul!(M, expτAQ, expτCQꜜ)
                for ϵ in ϵrange
                    # update q
                    updateq!(q, binomialcoefficients, uniquecoverages, n, ϵ)
                    loglik = 0.0
                    for index in eachindex(reads, counts)
                        loglik += calcloglik(counts, q, reads, coverages, M, hvals, index)
                    end
                    ns[index] = n 
                    τCs[index] = τC 
                    τAs[index] = τA
                    ϵs[index] = ϵ
                    logliks[index] = loglik 
                    if messages > 0 # No messages at all when messages ≤ 0.
                        if index % messages == 0
                            println("I did ", index, " posterior calculations, ", round(Int, 100*index / nrows, RoundDown), "% finished!")
                        end
                    end
                    index += 1 
                end # ϵ
            end # τA
        end # τC 
    end # n
    @assert index == nrows + 1 
    return (ns, τCs, τAs, ϵs, logliks)
end # function.
##
##
## End fixed R with counts
##

function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, R::Union{Integer, AbstractVector{<:Integer}}, reads::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=length(n)*length(τCrange)*length(τArange)*length(ϵrange)÷100)
    # Check conditions.
    all(0 .≤ ϵrange .<0.5) || throw(DomainError("ϵ should be in [0.0,0.5)."))
    length(nrange) > 0 || throw(ArgumentError("nrange should have at least one element."))
    minimum(nrange) > 0 || throw(ArgumentError("all n values in nrange should be positive."))
    length(τCrange) > 0 || throw(ArgumentError("τCrange should have at least one element."))
    minimum(τCrange) > 0 || throw(ArgumentError("All values in τCrange should be positive."))
    length(τArange) > 0 || throw(ArgumentError("τArange should have at least one element."))
    minimum(τArange) > 0 || throw(ArgumentError("All values in τArange should be positive."))
    length(ϵrange) > 0 || throw(ArgumentError("ϵrange should have at least one element."))
    minimum(ϵrange) > 0 || throw(ArgumentError("All values in ϵrange should be positive."))
    maximum(R) > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    minimum(reads) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(reads .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    checklength(y, reads, R) || throw(DimensionMismatch("y, reads (and, in case R is a vector, also R) should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0 and zero coverages. 
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with positive coverage should be strictly between 0 and 1."))
    
    y, reads, R, counts = filtervectorsandapplycountmap(y, reads, R, allowedindices)

    uniqueR = myunique(R)

    return exactposterior(nrange, τCrange, τArange, ϵrange, reads, R, uniqueR, counts, y, messages)
end # function.

#=
DICE 2-pop format:
ancestral_reads derived_reads frequency counts 

Can easily be loaded with 
CSV.read("filename.txt", DataFrame)
coveragevec = df[!, 1] + df[!, 2]
derivedreadsvec = df[!, 2]
freqs = df[!, 3]
counts = df[!, 4]
=#
"""
 MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, df::AbstractDataFrame, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

df is a DataFrame (from the DataFrames package) in DICE 2-pop format: It should contain exactly 4 columns, the first one contains the number of ancestral reads, the second the number of derived reads, the third the frequency of the derived allele in the anchor population and the fourth column contains the number of loci with exactly this combination of data. 
"""
function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, df::AbstractDataFrame, messages::Integer=nsteps÷100, scalingmessages::Bool=true)
    R = df[!, 1] + df[!, 2]
    d = df[!, 2]
    y = df[!, 3]
    counts = df[!, 4]
    return MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, R, d, y, counts, messages, scalingmessages)
end 

"""
 MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::AbstractString, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

dicefile is a (zipped) CSV-file in DICE 2-pop format: It should contain exactly 4 columns, the first one contains the number of ancestral reads, the second the number of derived reads, the third the frequency of the derived allele in the anchor population and the fourth column contains the number of loci with exactly this combination of data.
"""
function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::AbstractString, messages::Integer=nsteps÷100, scalingmessages::Bool=true)
    # Can handle both CSV file with and without headers. 
    df = readcsvfile(dicefile)
    return MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, df, messages, scalingmessages)
end  

"""
 MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::IO, messages::Integer=nsteps÷100, scalingmessages::Bool=true, header::Bool=true)

dicefile is an opened CSV-file in DICE 2-pop format: It should contain exactly 4 columns, the first one contains the number of ancestral reads, the second the number of derived reads, the third the frequency of the derived allele in the anchor population and the fourth column contains the number of loci with exactly this combination of data.

With the named argument header one can indicate whether the CSV file has a header. 
"""
function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::IO, messages::Integer=nsteps÷100, scalingmessages::Bool=true, header::Bool=true)
    # Can handle both CSV file with and without headers. 
    df = CSV.read(dicefile, DataFrame, header = header) 
    return MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, df, messages, scalingmessages)
end  


function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Union{Integer, AbstractVector{<:Integer}}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)
    checklength(y, d, R) || throw(DimensionMismatch("d, y  and/or R (in case R is a vector) should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with corresponding positive coverage should be strictly between 0 and 1."))

    y, d, R, counts = filtervectorsandapplycountmap(y, d, R, allowedindices)
    
    return MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, R, d, y, counts, messages, scalingmessages)
end 

"""
 MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Union{Integer, AbstractVector{<:Integer}}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, counts::AbstractVector{<:Integer}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

MCMC with simulated proposals. 
"""
function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::AbstractVector{<:Integer}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, counts::AbstractVector{<:Integer}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("For the current implementation we require that the prior on n is bounded by 10. I.e. P(prioronn ≤ 10)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but $prioronϵ can take negative values."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but $prioronϵ can take values larger than 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    any(r -> r>0, R) || throw(DomainError("At least some coverages should be positive.")) 
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    checklength(y, d, R, counts) || throw(DimensionMismatch("d, y, counts  and/or R (in case R is a vector) should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with corresponding positive coverage should be strictly between 0 and 1."))

    y, d, R, counts = filtervectorsandapplycountmap(y, d, R, counts, allowedindices)
    uniqueR = myunique(R)
    smallestn, largestn = extrema(prioronn)
    smallestτC, smallestτA = minimum(prioronτCτA)
    largestτC, largestτA = maximum(prioronτCτA)
    smallestϵ, largestϵ = extrema(prioronϵ)
    gridτC = range(smallestτC, largestτC, length=23)[begin+1:end] # avoid a 0.0 value
    gridτA = range(smallestτA, largestτA, length=23)[begin+1:end] # avoid a 0.0 value
    gridϵ = range(smallestϵ, largestϵ, length=24)[begin+1:end-1] # avoid a 0.0 and 0.5 value
    # ook nog aanpassen. 
    ns, τCs, τAs, ϵs, logliks = exactposterior(smallestn:largestn, gridτC, gridτA, gridϵ, d, R, uniqueR, counts, y, 0)
    maxlikindex = argmax(logliks)
    nmaxlik = ns[maxlikindex]
    τAmaxliksconditionedonn = Vector{Float64}(undef, nmax)
    τCmaxliksconditionedonn = Vector{Float64}(undef, nmax)
    ϵmaxliksconditionedonn = Vector{Float64}(undef, nmax)
    for n in 1:nmax 
        filteronn = ns .== n 
        τCsconditionedonn = τCs[filteronn]
        τAsconditionedonn = τAs[filteronn]
        ϵsconditionedonn = ϵs[filteronn]
        logliksconditionedonn = logliks[filteronn]
        maxlikindex = argmax(logliksconditionedonn)
        τCmaxliksconditionedonn[n] = τCsconditionedonn[maxlikindex]
        τAmaxliksconditionedonn[n] = τAsconditionedonn[maxlikindex]
        ϵmaxliksconditionedonn[n] = ϵsconditionedonn[maxlikindex]
        @assert unique(ns[filteronn]) == [n]
        @assert n ∉ ns[.!filteronn]
    end 
    
    d, R = makesizedvectors(d, R)
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax) # contains the matrices with values y^k(1-y)^{2n-k}.
    binomialcoefficients = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n] = preparedata(n, y)
    end
    hvals = [item for item in hvals] # change type of hvals to more specific one then Any.
    binomialcoefficients = [item for item in binomialcoefficients] # same.
    # Each item in results contains a tuple 
    # (nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob), 
    # which contain the states of the Markov chain. 
    results = Vector{Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, BitMatrix, Matrix{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        #p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = makeq(nmax, d, R)
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Matrix{Float64}(undef, nsteps, nmax)
        τAsample = Matrix{Float64}(undef, nsteps, nmax)
        ϵsample = Matrix{Float64}(undef, nsteps, nmax)
        acceptedvec = BitMatrix(undef, nsteps, nmax+1) # BitMatrix is more efficient than Matrix{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Matrix{Float64}(undef, nsteps, nmax+1)
        # Initialise.
        nprevious = nmaxlik
        τCprevious = copy(τCmaxliksconditionedonn)
        τAprevious = copy(τAmaxliksconditionedonn)
        ϵprevious = copy(ϵmaxliksconditionedonn)
        
        nsample[1] = nprevious
        τCsample[1,:] = τCprevious
        τAsample[1,:] = τAprevious 
        ϵsample[1,:] = ϵprevious
        acceptedvec[1,:] .= true 
        # Calculate log likelihood and log prior density with initial values.
        loglikprevious = Vector{Float64}(undef, nmax)
        loglikproposal = Vector{Float64}(undef, nmax)
        for n in 1:nmax
            loglikprevious[n] = logprobderivedreads!(q[n], d, R, uniqueR, counts, τCprevious[n], τAprevious[n], ϵprevious[n], n, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
        end
        # check if the initial value is valid. 
        any(loglikprevious .== -Inf) && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = Vector{Float64}(undef, nmax)
        for n in 1:nmax
            logpdfpriorprevious[n] = logpdf(prioronτCτA, [τCprevious[n], τAprevious[n]]) + logpdf(prioronϵ, ϵprevious[n])
        end
        logpdfnprev = logpdf(prioronn, nprevious)
                
        any(logpdfpriorprevious .== -Inf) && throw(DomainError("The initial value is outside the support of the prior."))
        logpdfnprev == -Inf && throw(DomainError("The initial value is outside the support of the prior on n."))
        logjointprob[1,2:end] .= loglikprevious .+ logpdfpriorprevious
        logjointprob[1,1] = loglikprevious[nprevious] + logpdfpriorprevious[nprevious] + logpdfnprev
        # Define proposal densities. 
        Σ = Diagonal(fill(0.000001, 2))
        distproposalτCτA = fill(MvNormal(Σ), nmax)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.0001
        distproposalϵ = fill(Normal(0.0, σ),nmax)
        nproposalprobvec = fill(0.2/(nmax-1), nmax)
        nproposalprobvec[nprevious] = 0.8
        distproposaln = Categorical(nproposalprobvec) # so with 80% chance nproposal is 
        scalefactor = fill(2.38^2/4, nmax) # = 2.38^2/4, as we have dimension of prior equal to 4.
        τCproposal = Vector{Float64}(undef, nmax)
        τAproposal = Vector{Float64}(undef, nmax)
        ϵproposal = Vector{Float64}(undef, nmax)
        logpdfpriorproposal = Vector{Float64}(undef, nmax)
        logacceptanceprobability = Vector{Float64}(undef, nmax)
        acceptanceprobability = Vector{Float64}(undef, nmax)
        accepted = BitVector(undef, nmax)
        for step in 2:nsteps  
            # draw proposals.
            for index in eachindex(τCproposal, τAproposal, ϵproposal) 
                τCproposal[index], τAproposal[index] = [τCprevious[index], τAprevious[index]] + rand(distproposalτCτA[index])
                ϵproposal[index] = ϵprevious[index] + rand(distproposalϵ[index])
            end
            for n in 1:nmax
                loglikproposal[n] = logprobderivedreads!(q[n], d, R, uniqueR, counts, τCproposal[n], τAproposal[n], ϵproposal[n], n, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
            end 
            for n in 1:nmax
                logpdfpriorproposal[n] = logpdf(prioronτCτA, [τCproposal[n], τAproposal[n]]) + logpdf(prioronϵ, ϵproposal[n])
            end 
            logacceptanceprobability .= loglikproposal .- loglikprevious .+ logpdfpriorproposal .- logpdfpriorprevious
            acceptanceprobability .= exp.(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            for index in eachindex(accepted, acceptanceprobability)
                accepted[index] = rand() ≤ acceptanceprobability[index]
            end 
            for n in 1:nmax
                if accepted[n]
                    τCsample[step, n] = τCprevious[n] = τCproposal[n]
                    τAsample[step, n] = τAprevious[n] = τAproposal[n]
                    ϵsample[step, n] = ϵprevious[n] = ϵproposal[n]
                    acceptedvec[step, n+1] = true
                    logjointprob[step, n+1] = loglikproposal[n] + logpdfpriorproposal[n]
                    # Update everything for the next cycle. 
                    loglikprevious[n] = loglikproposal[n]
                    logpdfpriorprevious[n] = logpdfpriorproposal[n]
                else
                    τCsample[step, n] = τCprevious[n]
                    τAsample[step, n] = τAprevious[n]
                    ϵsample[step, n] = ϵprevious[n]
                    acceptedvec[step, n+1] = false
                    logjointprob[step, n+1] = loglikprevious[n] + logpdfpriorprevious[n]
                    # no updating: loglikprevious and logpdfpriorprevious remain the same.
                end 
            end
            nproposal = rand(distproposaln) # so with 80% chance nproposal is equal to nprevious
            if nproposal == nprevious
                nsample[step] = nproposal
                acceptedvec[step, 1] = true
                logjointprob[step, 1] = logjointprob[step, nproposal+1] + logpdfnprev
            else 
                logpdfnprop = logpdf(prioronn, nproposal)
                logacceptprobn = logjointprob[step, nproposal+1] - logjointprob[step, nprevious+1] + logpdfnprop - logpdfnprev
                acceptprobn = exp(logacceptprobn)
                if rand() ≤ acceptprobn # proposal accepted 
                    nsample[step] = nprevious = nproposal
                    acceptedvec[step, 1] = true
                    logjointprob[step, 1] = logjointprob[step, nproposal+1] + logpdfnprop
                    nproposalprobvec = fill(0.2/(nmax-1), nmax)
                    nproposalprobvec[nprevious] = 0.8
                    distproposaln = Categorical(nproposalprobvec) # update the proposal distribution for n.
                    logpdfnprev = logpdfnprop
                else # proposal rejected 
                    nsample[step] = nprevious
                    acceptedvec[step, 1] = false 
                    logjointprob[step, 1] = logjointprob[step, nprevious+1] + logpdfnprev
                end 
            end                         
            if step % 1_000 == 0 && step ≥ 5_000
                # the covariance is calculated from the previous steps and used in the next iteration.
                for n in 1:nmax
                    Σ = Diagonal([var(τCsample[3_000:20:step, n], corrected=true), var(τAsample[3_000:20:step, n], corrected=true)])
                    σ = std(ϵsample[3_000:20:step, n], corrected=true)
                    if 10_000 < step < nsteps÷4  
                        acceptancerate = mean(acceptedvec[step-999:step, n+1]) # acceptance rate in the last 1000 steps.
                        if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                            scalefactor[n] = 3*scalefactor[n]/4
                            scalingmessages && println("Chain ", chainindex, " talking, I scaled the scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                        elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                            scalefactor[n] = 4*scalefactor[n]/3
                            scalingmessages && println("Chain ", chainindex, " talking, I scaled the scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                        end # So do nothing when the acceptancerate is between 0.1 and 0.6.     
                    end 
                    # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                    #Σ = scalefactor[n] * (cov([τCsample[1:10:step, n] τAsample[1:10:step, n]], corrected=true) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                    # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                    distproposalτCτA[n] = MvNormal(scalefactor[n]*Σ+1.0e-14*I)
                    # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                    # so not really necessary. 
                    distproposalϵ[n] = Normal(0.0, sqrt(scalefactor[n]) * σ + 1.0e-14)
                end 
            end
            if messages > 0 # No messages at all when messages ≤ 0.
                if step % messages == 0
                    println("Chain ", chainindex, " talking, I finished step ", step, " of ", nsteps, ", ", round(Int, 100*step / nsteps, RoundDown), " percent finished.")
                end
            end
        end # MCMC step
        results[chainindex] = (nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob)
    end # Chains. 
    return results
end # function.

end # module 

