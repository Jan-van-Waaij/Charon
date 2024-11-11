"""
    module PackageToBeShipped

This package provides functions to sample from the posterior. 
"""
module PackageToBeShipped

import Base: ==, exp, show # redefine exp and == for EigenExpansion types. 
# This requires an import statement. Redefine show for customized printing of EigenExpansion types. 

# external packages that I use. 
using Distributions # for probability distributions
using StatsBase # for the countmap function
using LinearAlgebra # Needed for the identity matrix I and the different matrix types and functions. 
using StaticArrays # for SizedVector, they are generally faster 
#than the build-in Matrix types, because they store their size in their type.
using Base.Threads # for multithreading.
using DataFrames # for processing Dice files. 
using CSV # for processing CSV files. 
      
# functions that you can use when you type 'using PackageToBeShipped'. 
export MCMCsampler # MCMC sampler.
export exactposterior # Calculate the posterior up to a fixed constant depending on the data, 
# but not depending on the parameters. Can be used for maximum likelihood estimation.
export unpackposterior # Constructs the posterior from the posterior from from the conditional chains. 

# Likelihood functions
# Tests written
"""
    EigenExpansion{Q<:AbstractMatrix{<:Real}, R<:Real, S<:AbstractVector{R}, T<:AbstractMatrix{<:Real}}

Representation of a matrix decomposition M=P*D*P^{-1}, where D is a diagonal matrix and P an invertible matrix with inverse Pinv. The fields are P, D and Pinv of types Q, Diagonal{R, S}, and T, respectively. 

This decomposition can be used for fast matrix exponentiation, using that exp(M) = P*exp(D)*Pinv, and the fact that the exponent of a diagonal matrix is formed by exponentiating its diagonal elements.
"""
struct EigenExpansion{Q<:AbstractMatrix{<:Real}, R<:Real, S<:AbstractVector{R}, T<:AbstractMatrix{<:Real}}
    P::Q
    D::Diagonal{R, S}
    Pinv::T
    # inner constructor.
    function EigenExpansion(M::AbstractMatrix{<:Real})
        F = eigen(M)
        P = F.vectors
        D = Diagonal(F.values)
        Pinv = inv(P)
        return new{typeof(P), eltype(F.values), typeof(F.values), typeof(Pinv)}(P, D, Pinv)
    end
end 

# outer constructor
# Tests written
"""
    EigenExpansion(M::Tridiagonal{S, T}}) where {S<:Real, T<:AbstractVector{S}}
    EigenExpansion(M::Diagonal{S, T}}) where {S<:Real, T<:AbstractVector{S}}
    EigenExpansion(M::SymTridiagonal{S, T}}) where {S<:Real, T<:AbstractVector{S}}
    EigenExpansion(M::Symmetric{S, T}}) where {S<:Real, T<:AbstractVector{S}}
    EigenExpansion(M::AbstractMatrix{<:Real})
    EigenExpansion(M::Bidiagonal{S, T}}) where {S<:Real, T<:AbstractVector{S}}

Try to find a decomposition M = P*D*P^(-1).
"""
EigenExpansion(M::Tridiagonal{S, T}) where {S<:Real, T<:AbstractVector{S}} = EigenExpansion(SizedMatrix{size(M)...,S}(M))
EigenExpansion(M::Diagonal{S, T}) where {S<:Real, T<:AbstractVector{S}} = EigenExpansion(SizedMatrix{size(M)...,S}(M))
EigenExpansion(M::Bidiagonal{S, T}) where {S<:Real, T<:AbstractVector{S}} = EigenExpansion(SizedMatrix{size(M)...,S}(M))
EigenExpansion(M::SymTridiagonal{S, T}) where {S<:Real, T<:AbstractVector{S}} = EigenExpansion(SizedMatrix{size(M)...,S}(M))
EigenExpansion(M::Symmetric{S, T}) where {S<:Real, T<:AbstractVector{S}} = EigenExpansion(SizedMatrix{size(M)...,S}(M))


"""
    Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::EigenExpansion{Q, R, S, T}) where {Q<:AbstractMatrix{<:Real}, R<:Real, S<:AbstractVector{R}, T<:AbstractMatrix{<:Real}}

Pretty display of EigenExpansion types. 
"""
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::EigenExpansion{Q, R, S, T}) where {Q<:AbstractMatrix{<:Real}, R<:Real, S<:AbstractVector{R}, T<:AbstractMatrix{<:Real}}
    print(io, "EigenExpansion{$Q, $R, $S, $T}\n\nP:\n\n")
    show(io, mime, F.P) 
    print(io, "\n\nD:\n\n")
    show(io, mime, F.D)
    print(io, "\n\nPinv:\n\n")
    show(io, mime, F.Pinv)
end



# Tests written. 
"""
    ==(E1::EigenExpansion, E2::EigenExpansion)

Test whether E1 and E2 have the same decomposition P*D*P^(-1), so whether E1.P == E2.P, E1.D == E2.D and E1.Pinv == E2.Pinv.
"""
==(E1::EigenExpansion, E2::EigenExpansion) = E1.P == E2.P && E1.D == E2.D && E1.Pinv == E2.Pinv

# tests written
"""
    exp(F::EigenExpansion, s::Real)

Calculate e^(s*M) in a numerical efficient way. 

If F=EigenExpansion(M), then exp(s*M)=exp(F,s)=F.P*exp(s*F.D)*F.Pinv, which makes 
use of the fast implementation of exp for diagonal matrices (it is just the exponetiation of the diagonal elements).     
"""
exp(F::EigenExpansion, s::Real) =  F.P*exp(s*F.D)*F.Pinv

# Tests written.
"""
    makeQ(n::Integer)

Calculate the (2n+1)x(2n+1) Q tridiagonal matrix as in Schraiber 2018. 
"""
function makeQ(n::Integer)
    lowerdiagonal = 0.5*(1.0:2n).*(0.0:2n-1.0)
    diagonal = -(0.0:2n).*(2n:-1.0:0.0)
    upperdiagonal = 0.5*(2n:-1.0:1.0).*(2n-1:-1.0:0.0)
    return Tridiagonal(lowerdiagonal, diagonal, upperdiagonal)
end 

# Tests written.
"""
    makeQꜜ(n::Integer)

Calculate the (2n+1)x(2n+1) Qꜜ tridiagonal matrix as in Schraiber 2018.  
"""
function makeQꜜ(n::Integer)
    lowerdiagonal = 0.5*(1.0:2n).*(0.0:2n-1.0)
    diagonal = -(0.0:2n).*(2n+1.0:-1.0:1.0)
    upperdiagonal = 0.5*(2n+1:-1.0:2.0).*(2n:-1.0:1.0)
    return Tridiagonal(lowerdiagonal, diagonal, upperdiagonal)
end

# Tests written. 
"""
    calcmatrix(Q::EigenExpansion, Qꜜ::EigenExpansion, τC::Real, τA::Real)

Calculate the matrix `e^{Q*τ_A}e^{Qꜜτ_C}` as in Schraiber. Note that in general, for matrices, exp(A+B) ≠ exp(A)exp(B) (equality holds if and only if A and B commute, which is not the case for Q and Qꜜ). 
"""
calcmatrix(Q::EigenExpansion, Qꜜ::EigenExpansion, τC::Real, τA::Real) =  exp(Q, τA)*exp(Qꜜ, τC)

# Tests written 
"""
    preparedata(n::Integer, y::AbstractVector{<:Real})

Returns a named tuple (Q=Q, Qꜜ=Qꜜ, hvals=hvals, binomialcoefficients=binomialcoefficients), where 
Q and Qꜜ are EigenExpansions of Q and Qꜜ (as in Schraiber 2018), respectively. 
hvals is a (2n+1) x length(y) matrix of type SizedMatrix. SizedMatrix is used rather than SMatrix, because this matrix is so large, and 
SMatrix is slow for large matrices. The matrix hvals is defined as hvals[k+1,i] = y[i]^k*(1-y[i])^(2n-k) (note that Julia has 1-based indexing).
binomialcoefficients is a vector (of type SVector) containing the binomial coefficients 2n over k, k=0,...,2n.  
"""
function preparedata(n::Integer, y::AbstractVector{<:Real})
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
    logprobderivedreads!(q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, τC::Real, τA::Real, ϵ::Real, n::Integer, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{<:Integer}, hvals::AbstractMatrix{<:Real})

Calculate the log probabilty of the data, a.k.a. the likelihood. It mutates q, which is dictionary, and is used for intermediate calculations. 
"""
function logprobderivedreads!(q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, τC::Real, τA::Real, ϵ::Real, n::Integer, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{<:Integer}, hvals::AbstractMatrix{<:Real})
    if τC < 0 || τA < 0 || ϵ < 0 || ϵ > 1
        return -Inf
    else
        M = calcmatrix(Q, Qꜜ, τC, τA) 
        updateq!(q, binomialcoefficients, uniquecoverages, n, ϵ)
        loglik = 0.0
        for index in eachindex(reads, coverages, counts)
            loglik += calcloglik(counts, q, reads, coverages, M, hvals, index)
        end
        return loglik
    end
end
# End likelihood functions. 

# tests written 
"""
    filtervectorsandapplycountmap(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}, allowedindices::AbstractVector{Bool})

This function does the following:
y_filtered = y[allowedindices]
d_filtered = d[allowedindices]
R_filtered = R[allowedindices]

Next it counts how many times (y_filtered[i], d_filtered[i], R_filtered[i]) occurs, for all triples, and next returns four vectors (y_filtered_unique, d_filtered_unique, R_filtered_unique, counts), where each triple (y_filtered_unique[i], d_filtered_unique[i], R_filtered_unique[i]) is unique and occurs counts[i] times in zip(y_filtered, d_filtered, R_filtered).
"""
function filtervectorsandapplycountmap(y::AbstractVector{<:Real}, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}, allowedindices::AbstractVector{Bool})
    y = y[allowedindices]
    d = d[allowedindices]
    R = R[allowedindices]
    tuples = zip(d, R, y)
    cm = countmap(tuples)
    d = [key[1] for key in keys(cm)]
    R = [key[2] for key in keys(cm)]
    y = [key[3] for key in keys(cm)]
    counts = collect(values(cm))
    return y, d, R, counts
end 

# tests written 
"""
    makeqfixedn(n::Integer, uniquedandRs::AbstractVector{<:Tuple{Integer, Integer}})

Construct a dictionary with keys (d[i], R[i]) and values real vectors of length 2n+1.
"""
makeqfixedn(n::Integer, uniquedandRs::AbstractVector{<:Tuple{Integer, Integer}}) = Dict((d, R) => Vector{Float64}(undef, 2n+1) for (d, R) in uniquedandRs)

# test written.
"""
    makeqfixedn(n::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}) 
"""
function makeqfixedn(n::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    uniquedandRs = unique(zip(d, R))
    return makeqfixedn(n, uniquedandRs)
end 

# tests written 
"""
    makeq(nmax::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})

Construct a vector of length nmax, where the n-th item is a dictionary with keys (d[i], R[i]) and values real vectors of length 2n+1. 
"""
function makeq(nmax::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    uniquedandRs = unique(zip(d, R))
    return [makeqfixedn(n, uniquedandRs) for n in 1:nmax]
end 

# tests written 
"""
    updateq!(q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, binomialcoefficients::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, n::Integer, ϵ::Real)

Update the dictionary q, so that q[(d, R)][k+1] is the probability of d derived reads out of R when you have k derived alleles with n inidividuals. 
"""
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

"""
    calcloglik(counts::AbstractVector{<:Integer}, q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, M::AbstractMatrix{<:Real}, hvals::AbstractMatrix{<:Real}, index::Integer)

Calculate the log-likelihood at locus index. 
"""
calcloglik(counts::AbstractVector{<:Integer}, q::AbstractDict{<:Tuple{Integer, Integer}, <:AbstractVector{<:Real}}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, M::AbstractMatrix{<:Real}, hvals::AbstractMatrix{<:Real}, index::Integer) = counts[index]*log(max(0.0, dot(q[reads[index], coverages[index]], M, view(hvals, :, index))))

"""
    readcsvfile(filename::AbstractString)

Read CSV file and turns it into a DataFrame (from DataFrames package). It automatically detects whether the CSV file has a header or not. 
"""
function readcsvfile(filename::AbstractString)
    # read the first line to determine if it is a header or not. 
    dftest = CSV.read(filename, DataFrame, header = false, limit=1)  
    if typeof(dftest[1,1]) <: AbstractString 
        # the first item in the first row was a string, so the csv file has a header.
        df = CSV.read(filename, DataFrame, header = true)
    else
        # the first item in the first row was not a string, so the csv file has no header. 
        df = CSV.read(filename, DataFrame, header = false)
    end 
    return df 
end 

"""
    exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer)

Calculate the exact posterior up to a fixed unknown constant that depends on the data, but does not depend on the parameters. 
"""
function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer)
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

"""
    exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, R::AbstractVector{<:Integer}, reads::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=length(n)*length(τCrange)*length(τArange)*length(ϵrange)÷100)
"""
function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, R::AbstractVector{<:Integer}, reads::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=length(n)*length(τCrange)*length(τArange)*length(ϵrange)÷100)
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
    length(y) == length(reads) == length(R) || throw(DimensionMismatch("y, reads and R should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0 and zero coverages. 
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with positive coverage should be strictly between 0 and 1."))
    
    y, reads, R, counts = filtervectorsandapplycountmap(y, reads, R, allowedindices)

    uniqueR = unique(R)

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
    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::AbstractVector{<:Integer}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, counts::AbstractVector{<:Integer}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

MCMC with simulated proposals. There are counts[i] loci with coverage R[i], d[i] derived alleles and frequency y[i] in the anchor population. 
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
    length(y) == length(d) == length(R) == length(counts) || throw(DimensionMismatch("d, y, counts  and R should be of the same length."))
    
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with corresponding positive coverage should be strictly between 0 and 1."))

    y = y[allowedindices]
    d = d[allowedindices]
    R = R[allowedindices]
    counts = counts[allowedindices]

    uniqueR = unique(R)
    smallestn, largestn = extrema(prioronn)
    smallestτC, smallestτA = minimum(prioronτCτA)
    largestτC, largestτA = maximum(prioronτCτA)
    smallestϵ, largestϵ = extrema(prioronϵ)
    gridτC = range(smallestτC, largestτC, length=23)[begin+1:end] # avoid a 0.0 value
    gridτA = range(smallestτA, largestτA, length=23)[begin+1:end] # avoid a 0.0 value
    gridϵ = range(smallestϵ, largestϵ, length=24)[begin+1:end-1] # avoid a 0.0 and 0.5 value
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
    
    d = SizedVector{length(d), Int}(d)
    R = SizedVector{length(R), Int}(R)
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
        # equal to the previous value, or equal to one of the other values 
        # 1,..., nprevious-1, nprevious+1,...,nmax with probability 0.2/(nmax-1).  
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
            # accept or reject the proposal with probability max(1, acceptanceprobability), 
            # note that rand() samples from the uniform distribution.
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
                        end # So do nothing when the acceptancerate is between 0.1 and 0.15.     
                    end 
                    distproposalτCτA[n] = MvNormal(scalefactor[n]*Σ+1.0e-14*I)
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

"""
    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::AbstractVector{<:Integer}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)

R[i] is the coverage, d[i] the number of derived reads and y[i] the frequency of the derived allele in the anchor population at locus i. 
"""
function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::AbstractVector{<:Integer}, d::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Integer=nsteps÷100, scalingmessages::Bool=true)
    length(y) == length(d) == length(R) || throw(DimensionMismatch("d, y  and/or R (in case R is a vector) should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with corresponding positive coverage should be strictly between 0 and 1."))

    y, d, R, counts = filtervectorsandapplycountmap(y, d, R, allowedindices)
    
    return MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, R, d, y, counts, messages, scalingmessages)
end 

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
    MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::IO, messages::Integer=nsteps÷100, scalingmessages::Bool=true, header::Bool=true)

dicefile is an opened CSV-file in DICE 2-pop format: It should contain exactly 4 columns, the first one contains the number of ancestral reads, the second the number of derived reads, the third the frequency of the derived allele in the anchor population and the fourth column contains the number of loci with exactly this combination of data.

If the file is a zipped file, open it using a zip decompressor:
```julia
using CodecZlib
dicefile = GzipDecompressorStream(open("path/to/your/file.csv.gz"))
```
With the named argument "header" one should indicate whether the CSV file has a header. Default is true.
"""
function MCMCsampler(nchains::Integer, nsteps::Integer, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, dicefile::IO, messages::Integer=nsteps÷100, scalingmessages::Bool=true, header::Bool=true)
    df = CSV.read(dicefile, DataFrame, header = header) 
    return MCMCsampler(nchains, nsteps, prioronn, prioronτCτA, prioronϵ, df, messages, scalingmessages)
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


function unpackposterior(chains::AbstractVector{<:Tuple{AbstractVector{<:Integer}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}, AbstractMatrix{Bool}, AbstractMatrix{<:Real}}}, nsteps::Integer)
    results = DataFrame(nsample = Int[],
            τCsample = Float64[],
            τAsample = Float64[],
            ϵsample = Float64[],
            accepted = Bool[],
            logjointprob = Float64[],
            chainid = Int[]
            )
    for index in eachindex(chains)
        nsample, τCsample, τAsample, ϵsample, _, logjointprob = chains[index]
        τCsampleunpacked = Vector{Float64}(undef, nsteps)
        τAsampleunpacked = Vector{Float64}(undef, nsteps)
        ϵsampleunpacked = Vector{Float64}(undef, nsteps)
        for index in 1:nsteps
            τCsampleunpacked[index] = τCsample[index, nsample[index]]
            τAsampleunpacked[index] = τAsample[index, nsample[index]]
            ϵsampleunpacked[index] = ϵsample[index, nsample[index]]
        end 
        acceptedvecunpacked = BitVector(undef, nsteps)
        acceptedvecunpacked[1] = true 
        for index in 2:nsteps
            everythingthesame = (τCsampleunpacked[index] == τCsampleunpacked[index-1]) && (τAsampleunpacked[index] == τAsampleunpacked[index-1]) && (ϵsampleunpacked[index] == ϵsampleunpacked[index-1]) && (nsample[index] == nsample[index-1])
            acceptedvecunpacked[index] = !everythingthesame
        end 
        df = DataFrame(nsample = nsample,
                            τCsample = τCsampleunpacked,
                            τAsample = τAsampleunpacked,
                            ϵsample = ϵsampleunpacked,
                            accepted = acceptedvecunpacked,
                            logjointprob = logjointprob[:,1],
                            chainid = fill(Int(index), nsteps)
                            )
        results = vcat(results, df)
    end 
    return results
end 

end # module 

