"""
    module PosteriorSampling

This package provides functions to sample from the posterior. 
"""
module PosteriorSampling

# external packages that I use. 
using Distributions, # for probability distributions
      Likelihood, # for the likelihood
      StatsBase, # for countmap
      LinearAlgebra, # Needed for the identity matrix I. 
      StaticArrays, # for SizedVector and MVector, they are generally faster.
      Base.Threads # for multithreading.
      
# functions that you can use when you type 'using PosteriorSampling'. 
export MHsampler, # The mcmc sampler.
       MCMCSimProp, # More efficient mcmc sampler for low coverage with simulated proposal, at the cost of a longer run time. 
       exactposterior

function filtervectorsandapplycountmap(y::Vector{Float64}, d::Vector{Int}, R::Vector{Int}, allowedindices::BitVector)
    y, d, R = y[allowedindices], d[allowedindices], R[allowedindices]
    tuples = zip(d, R, y)
    cm = countmap(tuples)
    d = [key[1] for key in keys(cm)]
    R = [key[2] for key in keys(cm)]
    y = [key[3] for key in keys(cm)]
    counts = collect(values(cm))
    return y, d, R, counts
end 

function filtervectorsandapplycountmap(y::Vector{Float64}, d::Vector{Int}, R::Int, allowedindices::BitVector)
    y, d = y[allowedindices], d[allowedindices]
    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
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
    MMatrix{2n+1, R+1, Float64}(undef)
end 

function makeqfixedn(n::Integer, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    uniquedandRs = unique(zip(d, R))
    return Dict((d, R) => Vector{Float64}(undef, 2n+1) for (d, R) in uniquedandRs)
end 

function makeq(nmax::Int, d::AbstractVector{<:Integer}, R::Int)
    return [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
end 

function makeq(nmax::Int, d::AbstractVector{<:Integer}, R::AbstractVector{<:Integer})
    uniquedandRs = unique(zip(d, R))
    return [Dict((d, R) => Vector{Float64}(undef, 2n+1) for (d, R) in uniquedandRs) for n in 1:nmax]
end 

checklength(y::Vector{Float64}, d::Vector{Int}, R::Vector{Int}) = length(y) == length(d) == length(R)

checklength(y::Vector{Float64}, d::Vector{Int}, R::Int) = length(y) == length(d)

function updateq!(q::AbstractMatrix{<:Real}, binomialcoefficients::AbstractVector{<:Integer}, uniquecoverages::Integer, n::Integer, ϵ::Real)
    for k in 0:2n
        dist = Binomial(uniquecoverages, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
        for d in 0:uniquecoverages
            q[k+1, d+1] = binomialcoefficients[k+1] * pdf(dist, d)
        end 
    end              
end 

function updateq!(q::AbstractDict, binomialcoefficients::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, n::Integer, ϵ::Real)
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

calcloglik(counts, q::AbstractMatrix, reads, coverages, M, hvals, index) = counts[index] * log(max(0.0, dot(view(q, :, reads[index]+1), M, view(hvals, :, index))))

calcloglik(counts, q::AbstractDict, reads, coverages, M, hvals, index) = counts[index]*log(max(0.0, dot(q[reads[index], coverages[index]], M, view(hvals, :, index))))

myunique(R::Integer) = R 

myunique(R::AbstractVector{<:Integer}) = unique(R)

##
##
## fixed R, error model, fixed n. 
##
##
"""
    exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::Union{AbstractVector{<:Integer}, Integer}, uniquecoverages::Union{AbstractVector{<:Integer}, Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Int)

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
function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::Union{AbstractVector{<:Integer}, Integer}, uniquecoverages::Union{AbstractVector{<:Integer}, Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Int)
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

function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, R::Union{Integer, AbstractVector{<:Integer}}, reads::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Int=length(n)*length(τCrange)*length(τArange)*length(ϵrange)÷100)
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

#
#
# begin variable R. 
#
#
#=
function exactposterior(nrange::AbstractVector{<:Integer}, τCrange::AbstractVector{<:Real}, τArange::AbstractVector{<:Real}, ϵrange::AbstractVector{<:Real}, reads::AbstractVector{<:Integer}, coverages::AbstractVector{<:Integer}, uniquecoverages::AbstractVector{<:Integer}, counts::AbstractVector{<:Integer}, y::AbstractVector{<:Real}, messages::Int)
    # Check conditions.
    all(0 .≤ ϵrange .< 0.5) || throw(DomainError("ϵ should be in [0.0,0.5)."))
    length(nrange) > 0 || throw(ArgumentError("nrange should have at least one element."))
    length(τCrange) > 0 || throw(ArgumentError("τCrange should have at least one element."))
    length(τArange) > 0 || throw(ArgumentError("τArange should have at least one element."))
    length(ϵrange) > 0 || throw(ArgumentError("ϵrange should have at least one element."))
    
    nrows = length(nrange)*length(τCrange)*length(τArange)*length(ϵrange)
    results = Matrix{Float64}(undef, 5, nrows)
    columnindex = 1

    for n in nrange
        Q, Qꜜ, hvals, binomialcoefficients = preparedata(n, y)
        #p = MVector{2n+1, Float64}(undef)
        q = makeqfixedn(n, reads, coverages)
        M = Matrix{Float64}(undef, 2n+1, 2n+1)
        for τC in τCrange 
            expτCQꜜ = exp(Qꜜ, τC)
            for τA in τArange
                expτAQ = exp(Q, τA)    
                mul!(M, expτAQ, expτCQꜜ)
                for ϵ in ϵrange
                    # update q. 
                    updateq!(q, uniquecoverages, n, ϵ)
                    loglik = 0.0
                    for index in eachindex(reads, coverages, counts)
                        # update p
                        #probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
                        #loglik +=  counts[index]*log(dot(p, q[reads[index], coverages[index]])) 
                        loglik += calcloglik(counts, q, reads, coverages, M, hvals, index)
                    end
                    results[:, columnindex] .= n, τC, τA, ϵ, loglik 
                    if messages > 0 # No messages at all when messages ≤ 0.
                        if columnindex % messages == 0
                            println("I did ", columnindex, " posterior calculations, ", round(Int, 100*columnindex / nrows, RoundDown), "% finished!")
                        end
                    end
                    columnindex += 1 
                end # ϵ
            end # τA
        end # τC 
    end # n
    @assert columnindex == nrows + 1 
    return results
end # function.
=#
#
#
# End variable R, with counts. 
#
#



##
##
## fixed R, error model, fixed n. 
##
##
"""
    MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing)

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
function MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100, scalingmessages::Bool=true)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0, 0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0, 0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    R > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    y = y[allowedindices]
    Rmax = R

    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
    
    smallestτC, smallestτA = minimum(prioronτCτA)
    largestτC, largestτA = maximum(prioronτCτA)
    smallestϵ, largestϵ = extrema(prioronϵ)
    gridτC = range(smallestτC, largestτC, length=23)[begin+1:end] # avoid a 0.0 value
    gridτA = range(smallestτA, largestτA, length=23)[begin+1:end] # avoid a 0.0 value
    gridϵ = range(smallestϵ, largestϵ, length=24)[begin+1:end-1] # avoid a 0.0 and 0.5 value

    _, τCs, τAs, ϵs, logliks = exactposterior([fixedn], gridτC, gridτA, gridϵ, d, R, uniqueR, counts, y, 0)

    maxlikindex = argmax(logliks)
    τCmaxlik = τCs[maxlikindex]
    τAmaxlik = τAs[maxlikindex]
    ϵmaxlik = ϵs[maxlikindex]
    
    d = SizedVector{length(d), Int}(d)
    counts = SizedVector{length(counts), Int}(counts)
    
    Q, Qꜜ, hvals, binomialcoefficients = preparedata(fixedn, y)
    results = Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, BitVector, Vector{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = MVector{2fixedn+1, Float64}(undef)
        q = MMatrix{2fixedn+1, R+1, Float64}(undef)
        #nsample = Vector{Int}(undef, nsteps)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        ϵsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps) # BitVector is more efficient than Vector{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Vector{Float64}(undef, nsteps)
        
        # Initialise.
        τCprevious = τCmaxlik 
        τAprevious = τAmaxlik
        ϵprevious = ϵmaxlik
        
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true
        loglikprevious = logprobderivedreads!(p, q, d, R, Rmax, counts, τCprevious, τAprevious, ϵprevious, fixedn, Q, Qꜜ, binomialcoefficients, hvals)
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = loglikprevious + logpdfpriorprevious
        # Define proposal densities. 
        Σ = Diagonal(fill(0.000001, 2))
        distproposalτCτA = MvNormal(Σ)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.0001
        distproposalϵ = Normal(0.0, σ)
        scalefactor = 2.38^2/3.0 #0.1 #1.8881 # = 2.38^2/3, as we have dimension of prior equal to 3.
        for step in 2:nsteps  
            # draw proposals. 
            τCproposal, τAproposal = [τCprevious, τAprevious] + rand(distproposalτCτA)
            ϵproposal = ϵprevious + rand(distproposalϵ)
            # calculate the loglikelihood
            loglikproposal = logprobderivedreads!(p, q, d, R, Rmax, counts, τCproposal, τAproposal, ϵproposal, fixedn, Q, Qꜜ, binomialcoefficients, hvals)
            logpdfpriorproposal = logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = loglikproposal - loglikprevious + logpdfpriorproposal - logpdfpriorprevious 
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = loglikproposal + logpdfpriorproposal
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = loglikprevious + logpdfpriorprevious 
                # no updating: loglikprevious and logpdfpriorprevious remain the same.
            end 
            if step % 1000 == 0 && step ≥ 5000
                # the covariance is calculated from the previous steps and used in the next iteration.
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = Diagonal([var(τCsample[3000:20:step], corrected=true), var(τAsample[3000:20:step], corrected=true)])# + 1.0e-14*I
                σ = std(ϵsample[3000:20:step], corrected=true) # + 1.0e-14     
                # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                if 10_000 < step < nsteps÷4   
                    acceptancerate = mean(acceptedvec[step-999:step]) # acceptance rate in the last 1000 steps.
                    if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                        scalefactor = 3*scalefactor/4
                        scalingmessages && println("Chain ", chainindex, " talking, I scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                    elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                        scalefactor = 4*scalefactor/3
                        scalingmessages && println("Chain ", chainindex, " talking, I scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                    end # So do nothing when the acceptancerate is between 0.1 and 0.6. 
                end
                distproposalτCτA = MvNormal(scalefactor*Σ+1.0e-14*I)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                # small threshold to prevent getting zero variance. 
                distproposalϵ = Normal(0.0, sqrt(scalefactor)*σ + 1.0e-14) 
            end
            if messages > 0 # No messages at all when messages ≤ 0.
                if step % messages == 0
                    println("Chain ", chainindex, " talking, I finished step ", step, " of ", nsteps, ", ", round(Int, 100*step / nsteps, RoundDown), " percent finished.")
                end
            end
        end # MCMC step
        results[chainindex] = (τCsample, τAsample, ϵsample, acceptedvec, logjointprob)
    end # Chains. 
    return results
end # function.
##
##
## End fixed R, error model, fixed n. 
##

# fixed R, error model, unknown n.
"""
    MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on n, τC, τA and ϵ. Instead of using the recommended scaling factor 1.4161, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress.

If not specified, the initial values are choosen from the prior.

This function returns a vector of length nchains. Each element contains the results of one chain; it is a tuple of six vectors: a vector containing the n-values, a vector containing the τC-values, a vector containing the τA-values, a vector containing the ϵ values, a vector telling whether a certain proposal was accepted or not and a vector with unnormalised posterior log-densities. So for example 

```julia
results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, R, d, y, messages)
nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[1]
```

gives the results of the first chain.

Here we assume that all coverages are the same.
"""
function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Union{Int, Vector{Int}}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100, scalingmessages::Bool=true)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("For the current implementation we require that the prior on n is bounded by 10. I.e. P(prioronn ≤ 10)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain, you provided $n chains."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    Rmax = maximum(R)
    Rmax > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    sum(allowedindices) > 0 || throw(ArgumentError("At least one frequency where the coverage is positive should be strictly between 0.0 and 1.0."))
    y, d, R, counts = filtervectorsandapplycountmap(y, d, R, allowedindices)
    uniqueR = myunique(R)

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
    τCmaxlik = τCs[maxlikindex]
    τAmaxlik = τAs[maxlikindex]
    ϵmaxlik = ϵs[maxlikindex]
    
    d, R = makesizedvectors(d, R)
    counts = SizedVector{length(counts), Int}(counts)

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
    results = Vector{Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}, Vector{Float64}, BitVector, Vector{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        #p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = makeq(nmax, d, R)
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        ϵsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps) # BitVector is more efficient than Vector{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Vector{Float64}(undef, nsteps)
        # Initialise.
        nprevious = nmaxlik
        τCprevious = τCmaxlik
        τAprevious = τAmaxlik
        ϵprevious = ϵmaxlik
        
        nsample[1] = nprevious
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious 
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true 
        # Calculate log likelihood and log prior density with initial values.
        # p niet meer nodig. 
        loglikprevious = logprobderivedreads!(q[nprevious], d, R, uniqueR, counts, τCprevious, τAprevious, ϵprevious, nprevious, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])

        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronn, nprevious) + 
            logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = loglikprevious + logpdfpriorprevious
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = MvNormal(Σ)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.001
        distproposalϵ = Normal(0.0, σ)
        scalefactor = 2.38^2/4 # = 2.38^2/4, as we have dimension of prior equal to 4.
        for step in 2:nsteps  
            # draw proposals. 
            τCproposal, τAproposal = [τCprevious, τAprevious] + rand(distproposalτCτA)
            # The distribution of the next line is symmetric, i.e. P(a;b)=P(b;a), 
            # even for a,b equal or near 1 or n.  
            nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            ϵproposal = ϵprevious + rand(distproposalϵ)
            # calculate the loglikelihood
            loglikproposal = logprobderivedreads!(q[nproposal], d, R, uniqueR, counts, τCproposal, τAproposal, ϵproposal, nproposal, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])

            logpdfpriorproposal = logpdf(prioronn, nproposal) + logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                nsample[step] = nprevious = nproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = loglikproposal + logpdfpriorproposal
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                nsample[step] = nprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = loglikprevious + logpdfpriorprevious 
                # no updating: loglikprevious and logpdfpriorprevious remain the same.
            end 
            if step % 1000 == 0 && step ≥ 5000
                Σ = Diagonal([var(τCsample[3000:20:step], corrected=true), var(τAsample[3000:20:step], corrected=true)])# + 1.0e-14*I
                σ = std(ϵsample[3000:20:step], corrected=true) # + 1.0e-14     
                # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                if 10_000 < step < nsteps÷4   
                    acceptancerate = mean(acceptedvec[step-999:step]) # acceptance rate in the last 1000 steps.
                    if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                        scalefactor = 3*scalefactor/4
                        scalingmessages && println("Chain ", chainindex, " talking, I scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                    elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                        scalefactor = 4*scalefactor/3
                        scalingmessages && println("Chain ", chainindex, " talking, I scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                    end # So do nothing when the acceptancerate is between 0.1 and 0.6. 
                end
                distproposalτCτA = MvNormal(scalefactor*Σ+1.0e-14*I)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                # small threshold to prevent getting zero variance. 
                distproposalϵ = Normal(0.0, sqrt(scalefactor)*σ + 1.0e-14) 
                # the covariance is calculated from the previous steps and used in the next iteration.
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
## end new, adding error model.

"""
    MCMCSimProp(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Union{Int, Vector{Int}}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100, scalingmessages::Bool=true)

MCMC with simulated proposals. 
"""
function MCMCSimProp(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Union{Int, Vector{Int}}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100, scalingmessages::Bool=true)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("For the current implementation we require that the prior on n is bounded by 10. I.e. P(prioronn ≤ 10)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    any(r -> r>0, R) || throw(DomainError("At least some coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    checklength(y, d, R) || throw(DimensionMismatch("d, y  and/or R (in case R is a vector) should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    any(allowedindices) || throw(ArgumentError("At least one element of y with corresponding positive coverage should be strictly between 0 and 1."))

    y, d, R, counts = filtervectorsandapplycountmap(y, d, R, allowedindices)
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

end # package PosteriorSampling
