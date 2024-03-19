"""
    module PosteriorSampling

This package provides functions to sample from the posterior. 
"""
module PosteriorSampling

using Distributions, # for probability distributions
      Likelihood, # for the likelihood
      StatsBase, # for countmap
      LinearAlgebra, # Needed for the identity matrix I. 
      StaticArrays, # for SizedVector and MVector, they are generally faster.
      Base.Threads # for multithreading.
      
export MHsampler, # The mcmc sampler.
       MCMCSimProp, # More efficient mcmc sampler for low coverage with simulated proposal, at the cost of a longer run time. 
       SpecialSampler, 
       ExactPosterior

##
##
## fixed R, error model, fixed n. 
##
##
"""
    ExactPosterior(fixedn::Int, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, fixedϵ::Float64, prioronτCτA::ContinuousMultivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100)

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
function ExactPosterior(fixedn::Int, τrange::AbstractVector{Float64}, fixedϵ::Float64, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=length(τrange)^2÷100)
    # Check conditions.
    0≤fixedϵ<0.5 || throw(DomainError("ϵ should be in [0.0,0.5)."))
    length(τrange) > 0 || throw(ArgumentError("τrange should have at least one element."))
    R > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    any(allowedindices) || throw(ArgumentError("At least one element of y should be strictly between 0 and 1."))
    d = d[allowedindices]
    y = y[allowedindices]

    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
    d = SizedVector{length(d), Int}(d)

    Q, Qꜜ, hvals, binomialcoefficients = preparedata(fixedn, y)
    p = MVector{2fixedn+1, Float64}(undef)
    q = MMatrix{2fixedn+1, R+1, Float64}(undef)
    M = Matrix{Float64}(undef, 2fixedn+1, 2fixedn+1)
        
    for k in 0:2n
        dist = Binomial(R, k*(1.0-fixedϵ)/(2n) + (2n-k)*fixedϵ/(2n))
        for d in 0:R
            q[k+1, d+1] = pdf(dist, d)
        end 
    end 
    
    nrows = length(τrange)^2
    results = Matrix{Float64}(undef, 3, nrows)
    columnindex = 1
    for τC in τrange 
        expτCQꜜ = exp(Qꜜ, τC)
        for τA in τrange
            expτAQ = exp(Q, τA)    
            mul!(M, expτAQ, expτCQꜜ)
            loglik = 0.0
            for index in eachindex(reads, counts)
                # update p
                probderivedalleles!(p, view(hvals, :, index), M, binomialcoefficients)
                loglik +=  counts[index]*log(dot(p, view(q, :, reads[index]+1))) 
            end
            results[:, columnindex] .= τC, τA, loglik 
            if messages > 0 # No messages at all when messages ≤ 0.
                if columnindex % messages == 0
                    println("I did ", columnindex, " posterior calculations, ", round(Int, 100*columnindex / nrows, RoundDown), "% finished!")
                end
            end
            columnindex += 1 
        end
    end 
    @assert columnindex == nrows + 1 
    return results
end # function.
##
##
## End fixed R, error model, fixed n. 
##

##
##
## fixed R, error model, fixed n. 
##
##
"""
    SpecialSampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing)

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
function SpecialSampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing, stdτ::Float64, stdϵ::Float64)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0, 0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0, 0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    #R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    R > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    #R = R[allowedindices]
    #R = SizedVector{length(R), Int}(R)
    y = y[allowedindices]

    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
    for index in eachindex(d, y, counts)
        @assert cm[(d[index], y[index])] == counts[index]
    end 
    cm2 = Dict((d[index], y[index]) => counts[index] for index in eachindex(d, y, counts))
    @assert cm == cm2 
    d = SizedVector{length(d), Int}(d)
    Q, Qꜜ, hvals, binomialcoefficients = preparedata(fixedn, y)

    Σ = stdτ^2*I(2)#10.0^(-10)*I(2)
    distproposalτCτA = MvNormal(Σ)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
    distproposalϵ = Normal(0.0, stdϵ)

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
        if init === nothing # draw initial value from the prior.
            #nprevious = rand(prioronn)
            τCprevious, τAprevious = rand(prioronτCτA)
            ϵprevious = rand(prioronϵ)
        else
            #nprevious, 
            τCprevious, τAprevious, ϵprevious = init 
        end
        #nsample[1] = nprevious
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true
        loglikprevious = logprobderivedreads!(p, q, d, counts, τCprevious, τAprevious, ϵprevious, fixedn, R, Q, Qꜜ, binomialcoefficients, hvals)
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = loglikprevious + logpdfpriorprevious 
        for step in 2:nsteps  
            # draw proposals. 
            τCproposal, τAproposal = [τCprevious, τAprevious] + rand(distproposalτCτA)
            # The distribution of the next line is symmetric, i.e. P(a;b)=P(b;a), 
            # even for a,b equal or near 1 or n.  
            #nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            #nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            ϵproposal = ϵprevious + rand(distproposalϵ)
            # calculate the loglikelihood
            loglikproposal = logprobderivedreads!(p, q, d, counts, τCproposal, τAproposal, ϵproposal, fixedn, R, Q, Qꜜ, binomialcoefficients, hvals)
            logpdfpriorproposal = logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
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
## End SpecialSampler 
##

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
function MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing, factor::Union{Nothing, Float64}=nothing, temperature::Float64=1.0)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0, 0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0, 0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    #R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    R > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    1.0 ≤ temperature < Inf || throw(DomainError("temperature should be larger than 1.0 and finite."))
    inversetemperature = 1/temperature

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    #R = R[allowedindices]
    #R = SizedVector{length(R), Int}(R)
    y = y[allowedindices]

    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
    for index in eachindex(d, y, counts)
        @assert cm[(d[index], y[index])] == counts[index]
    end 
    cm2 = Dict((d[index], y[index]) => counts[index] for index in eachindex(d, y, counts))
    @assert cm == cm2 
    d = SizedVector{length(d), Int}(d)



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
        if init === nothing # draw initial value from the prior.
            #nprevious = rand(prioronn)
            τCprevious, τAprevious = rand(prioronτCτA)
            ϵprevious = rand(prioronϵ)
        else
            #nprevious, 
            τCprevious, τAprevious, ϵprevious = init 
        end
        #nsample[1] = nprevious
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true
        loglikprevious = logprobderivedreads!(p, q, d, counts, τCprevious, τAprevious, ϵprevious, fixedn, R, Q, Qꜜ, binomialcoefficients, hvals)
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = loglikprevious + logpdfpriorprevious
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = MvNormal(Σ)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.001
        distproposalϵ = Normal(0.0, σ)
        if factor === nothing 
            scalefactor = 1.8881 # = 2.38^2/3, as we have dimension of prior equal to 3.
        else 
            scalefactor = factor * 1.8881
        end  
        for step in 2:nsteps  
            # draw proposals. 
            τCproposal, τAproposal = [τCprevious, τAprevious] + rand(distproposalτCτA)
            # The distribution of the next line is symmetric, i.e. P(a;b)=P(b;a), 
            # even for a,b equal or near 1 or n.  
            #nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            #nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            ϵproposal = ϵprevious + rand(distproposalϵ)
            # calculate the loglikelihood
            loglikproposal = logprobderivedreads!(p, q, d, counts, τCproposal, τAproposal, ϵproposal, fixedn, R, Q, Qꜜ, binomialcoefficients, hvals)
            logpdfpriorproposal = logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = inversetemperature * ( loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious )
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
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
                if factor === nothing 
                    if step < nsteps÷4 && step > 1000  
                        acceptancerate = mean(acceptedvec[step-999:step]) # acceptance rate in the last 1000 steps.
                        if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                            scalefactor = 3*scalefactor/4
                            println("Chain ", chainindex, " talking, I scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                        elseif acceptancerate > 0.6 # acceptancerate too large, take bigger steps.
                            scalefactor = 4*scalefactor/3
                            println("Chain ", chainindex, " talking, I scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                        end # So do nothing when the acceptancerate is between 0.1 and 0.6. 
                    end 
                end 
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]], corrected=true) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step], corrected=true) + 1.0e-10 # small threshold to prevent getting zero variance. 
                distproposalϵ = Normal(0.0, sqrt(scalefactor) * σ)
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
function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    maximum(prioronn) < Inf || throw(ArgumentError("For the current implementation we require that the prior on n is bounded. I.e. there is an N, so that P(prioronn ≤ N)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,1], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain, you provided $n chains."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    R > 0 || throw(DomainError("All coverages should be positive.")) 
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    maximum(d) ≤  R || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    sum(allowedindices) > 0 || throw(ArgumentError("At least one frequency should be strictly between 0.0 and 1.0."))
    d = d[allowedindices]
    y = y[allowedindices]
    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
    for index in eachindex(d, y, counts)
        @assert cm[(d[index], y[index])] == counts[index]
    end 
    cm2 = Dict((d[index], y[index]) => counts[index] for index in eachindex(d, y, counts))
    @assert cm == cm2 
    d = SizedVector{length(d), Int}(d)

    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("The current implementation can handle at most 10 individuals."))
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
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        ϵsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps) # BitVector is more efficient than Vector{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Vector{Float64}(undef, nsteps)
        # Initialise.
        if init === nothing # draw initial value from the prior.
            nprevious = rand(prioronn)
            τCprevious, τAprevious = rand(prioronτCτA)
            ϵprevious = rand(prioronϵ)
        else 
            nprevious, τCprevious, τAprevious, ϵprevious = init 
        end 
        nsample[1] = nprevious
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious 
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true 
        # Calculate log likelihood and log prior density with initial values.
        # HIER MOET IK DUS logprobderivedreads GAAN VERANDEREN!!!!
        loglikprevious = logprobderivedreads!(p[nprevious], q[nprevious], d, counts, τCprevious, τAprevious, ϵprevious, nprevious, R, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])
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
        scalefactor = 1.4161 # = 2.38^2/4, as we have dimension of prior equal to 4.
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
            loglikproposal = logprobderivedreads!(p[nproposal], q[nproposal], d, counts, τCproposal, τAproposal, ϵproposal, nproposal, R, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])
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
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
                if step < nsteps÷4 && step > 1000  
                    percentage = round(Int, 100*step / nsteps)
                    acceptancerate = mean(acceptedvec[step-999:step]) # acceptance rate in the last 1000 steps.
                    if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                        scalefactor = 3*scalefactor/4
                        println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                    elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                        scalefactor = 4*scalefactor/3
                        println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                    end # So do nothing when the acceptancerate is between 0.1 and 0.6. 
                end 
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]], corrected=true) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step], corrected=true) + 1.0e-5 # small threshold to prevent getting zero variance. 
                distproposalϵ = Normal(0.0, sqrt(scalefactor) * σ)
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

# variable R, error model, unknown n. 
"""
    MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on n, τC, τA and ϵ. Instead of using the recommended scaling factor 1.4161, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress.

If not specified, the initial values are choosen from the prior.

This function returns a vector of length nchains. Each element contains the results of one chain; it is a tuple of six vectors: a vector containing the n-values, a vector containing the τC-values, a vector containing the τA-values, a vector containing the ϵ values, a vector telling whether a certain proposal was accepted or not and a vector with unnormalised posterior log-densities. So for example 

```julia
results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, R, d, y, messages)
nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[1]
```

gives the results of the first chain.
"""
function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    maximum(prioronn) < Inf || throw(ArgumentError("For the current implementation we require that the prior on n is bounded. I.e. there is an N, so that P(prioronn ≤ N)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("The current implementation can handle at most 10 individuals."))
    Rmax = maximum(R)
    Rmax > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    d = d[allowedindices]
    R = R[allowedindices]
    y = y[allowedindices]
    tuples = zip(d, R, y)
    cm = countmap(tuples)
    d = [key[1] for key in keys(cm)]
    R = [key[2] for key in keys(cm)]
    y = [key[3] for key in keys(cm)]
    counts = collect(values(cm))
    for index in eachindex(d, R, y, counts)
        @assert cm[(d[index], R[index], y[index])] == counts[index]
    end 
    cm2 = Dict((d[index], R[index], y[index]) => counts[index] for index in eachindex(d, R, y, counts))
    @assert cm == cm2 

    d = SizedVector{length(d), Int}(d)
    R = SizedVector{length(R), Int}(R)
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax) # contains the matrices with values y^k(1-y)^{2n-k}.
    binomialcoefficients = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n] = preparedata(n, y)
    end
    #p = Vector{Any}(undef, nmax)
    #binomialprobs = Vector{Any}(undef, nmax)
    hvals = [item for item in hvals] # change type of hvals to more specific one then Any.
    binomialcoefficients = [item for item in binomialcoefficients] # same.
    # Each item in results contains a tuple 
    # (nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob), 
    # which contain the states of the Markov chain. 
    results = Vector{Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}, Vector{Float64}, BitVector, Vector{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = [[MMatrix{2n+1, R+1, Float64}(undef) for R in 1:Rmax] for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        ϵsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps) # BitVector is more efficient than Vector{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Vector{Float64}(undef, nsteps)
        # Initialise.
        if init === nothing # draw initial value from the prior.
            nprevious = rand(prioronn)
            τCprevious, τAprevious = rand(prioronτCτA)
            ϵprevious = rand(prioronϵ)
        else
            nprevious, τCprevious, τAprevious, ϵprevious = init 
        end
        nsample[1] = nprevious
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true
        loglikprevious = logprobderivedreads!(p[nprevious], q[nprevious], d, R, counts, Rmax, τCprevious, τAprevious, ϵprevious, nprevious, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])
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
        scalefactor = 1.4161 # = 2.38^2/4, as we have dimension of prior equal to 4.
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
            loglikproposal = logprobderivedreads!(p[nproposal], q[nproposal], d, R, counts, Rmax, τCproposal, τAproposal, ϵproposal, nproposal, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])
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
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
                if step < nsteps÷4 && step > 1000  
                    acceptancerate = mean(acceptedvec[step-999:step]) # acceptance rate in the last 1000 steps.
                    if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                        scalefactor = 3*scalefactor/4
                        println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                    elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                        scalefactor = 4*scalefactor/3
                        println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                    end # So do nothing when the acceptancerate is between 0.1 and 0.6. 
                end 
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]], corrected=true) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step], corrected=true) + 1.0e-5 # small threshold to prevent getting zero variance. 
                distproposalϵ = Normal(0.0, sqrt(scalefactor) * σ)
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
# End variable R, error model. 



##
##
## MCMC with simulated proposals. 
# function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
function MCMCSimProp(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    maximum(prioronn) < Inf || throw(ArgumentError("For the current implementation we require that the prior on n is bounded. I.e. there is an N, so that P(prioronn ≤ N)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    maximum(d) ≤  R || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("The current implementation can handle at most 10 individuals."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    y = y[allowedindices]

    tuples = zip(d, y)
    cm = countmap(tuples)
    thekeys = keys(cm)
    d = [key[1] for key in thekeys]
    y = [key[2] for key in thekeys]
    counts = collect(values(cm))
    for index in eachindex(d, y, counts)
        @assert cm[(d[index], y[index])] == counts[index]
    end 
    cm2 = Dict((d[index], y[index]) => counts[index] for index in eachindex(d, y, counts))
    @assert cm == cm2 
    
    d = SizedVector{length(d), Int}(d)
    #R = SizedVector{length(R), Int}(R)
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax) # contains the matrices with values y^k(1-y)^{2n-k}.
    binomialcoefficients = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n] = preparedata(n, y)
    end
    #p = Vector{Any}(undef, nmax)
    #binomialprobs = Vector{Any}(undef, nmax)
    hvals = [item for item in hvals] # change type of hvals to more specific one then Any.
    binomialcoefficients = [item for item in binomialcoefficients] # same.
    # Each item in results contains a tuple 
    # (nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob), 
    # which contain the states of the Markov chain. 
    results = Vector{Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, BitMatrix, Matrix{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Matrix{Float64}(undef, nsteps, nmax)
        τAsample = Matrix{Float64}(undef, nsteps, nmax)
        ϵsample = Matrix{Float64}(undef, nsteps, nmax)
        acceptedvec = BitMatrix(undef, nsteps, nmax+1) # BitMatrix is more efficient than Matrix{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Matrix{Float64}(undef, nsteps, nmax+1)
        # Initialise.
        τCprevious = Vector{Float64}(undef, nmax)
        τAprevious = Vector{Float64}(undef, nmax)
        ϵprevious = Vector{Float64}(undef, nmax)
        if init === nothing # draw initial value from the prior.
            nprevious = rand(prioronn)
            for index in eachindex(τCprevious, τAprevious)
                τCprevious[index], τAprevious[index] = rand(prioronτCτA)
                ϵprevious[index] = rand(prioronϵ)
            end 
        else 
            nprevious = init[1]
            for index in eachindex(τCprevious, τAprevious)
                _, τCprevious[index], τAprevious[index], ϵprevious[index] = init
            end
        end 
        nsample[1] = nprevious
        τCsample[1,:] = τCprevious
        τAsample[1,:] = τAprevious 
        ϵsample[1,:] = ϵprevious
        acceptedvec[1,:] .= true 
        # Calculate log likelihood and log prior density with initial values.
        loglikprevious = Vector{Float64}(undef, nmax)
        loglikproposal = Vector{Float64}(undef, nmax)
        for n in 1:nmax
            loglikprevious[n] = logprobderivedreads!(p[n], q[n], d, counts, τCprevious[n], τAprevious[n], ϵprevious[n], n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
        end
        # check if the initial value is valid. 
        any(loglikprevious .== -Inf) && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = Vector{Float64}(undef, nmax)
        for n in 1:nmax
            logpdfpriorprevious[n] = logpdf(prioronτCτA, [τCprevious[n], τAprevious[n]]) + logpdf(prioronϵ, ϵprevious[n])
        end
        any(logpdfpriorprevious .== -Inf) && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1,2:end] .= loglikprevious .+ logpdfpriorprevious
        logjointprob[1,1] = loglikprevious[nprevious] + logpdfpriorprevious[nprevious] + logpdf(prioronn, nprevious)
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = repeat([MvNormal(Σ)], nmax)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.001
        distproposalϵ = repeat([Normal(0.0, σ)],nmax)
        scalefactor = repeat([1.8881], nmax) # = 2.38^2/4, as we have dimension of prior equal to 4.
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
                loglikproposal[n] = logprobderivedreads!(p[n], q[n], d, counts, τCproposal[n], τAproposal[n], ϵproposal[n], n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
            end 
            for n in 1:nmax
                logpdfpriorproposal[n] = logpdf(prioronτCτA, [τCproposal[n], τAproposal[n]]) + logpdf(prioronϵ, ϵproposal[n])
            end 
            logacceptanceprobability .= loglikproposal .- loglikprevious .+ logpdfpriorproposal .- logpdfpriorprevious
            acceptanceprobability .= exp.(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            for index in eachindex(acceptanceprobability)
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
            nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            logpdfnprev = logpdf(prioronn, nprevious)
            if nproposal == nprevious
                nsample[step] = nproposal
                acceptedvec[step, 1] = true
                logjointprob[step, 1] = logjointprob[step, nproposal+1] + logpdfnprev
            else 
                logpdfnprop = logpdf(prioronn, nproposal)
                logacceptprobn = logjointprob[step, nproposal+1] - logjointprob[step, nprevious+1] + logpdfnprop - logpdfnprev
                acceptprobn = exp(logacceptprobn)
                if rand() ≤ acceptprobn
                    nsample[step] = nproposal
                    nprevious = nproposal
                    acceptedvec[step, 1] = true
                    logjointprob[step, 1] = logjointprob[step, nproposal+1] + logpdfnprop
                else
                    nsample[step] = nprevious
                    acceptedvec[step, 1] = false 
                    logjointprob[step, 1] = logjointprob[step, nprevious+1] + logpdfnprev
                end 
            end             
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
                for n in 1:nmax
                    if step < nsteps÷4 && step > 1000  
                        acceptancerate = mean(acceptedvec[step-999:step, n+1]) # acceptance rate in the last 1000 steps.
                        if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                            scalefactor = 3*scalefactor/4
                            println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                        elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                            scalefactor = 4*scalefactor/3
                            println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                        end # So do nothing when the acceptancerate is between 0.1 and 0.6.     
                    end 
                    # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                    Σ = scalefactor[n] * (cov([τCsample[1:10:step, n] τAsample[1:10:step, n]], corrected=true) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                    # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                    distproposalτCτA[n] = MvNormal(Σ)
                    # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                    # so not really necessary. 
                    σ = std(ϵsample[1:10:step, n], corrected=true) + 1.0e-10 # small threshold to prevent getting zero variance. 
                    distproposalϵ[n] = Normal(0.0, sqrt(scalefactor[n]) * σ)
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

##
##
## MCMCSimProp non-constant coverage. 
##
##
function MCMCSimProp(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    maximum(prioronn) < Inf || throw(ArgumentError("For the current implementation we require that the prior on n is bounded. I.e. there is an N, so that P(prioronn ≤ N)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the minimum is negative."))
    maximum(prioronϵ) ≤ 0.5 || throw(ArgumentError("prioronϵ should have support on [0,0.5], but the maximum exceeds 0.5."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    Rmax = maximum(R)
    Rmax > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) == length(R) || throw(DimensionMismatch("R, d and y should be of the same length."))
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    d = d[allowedindices]
    R = R[allowedindices]
    y = y[allowedindices]
    tuples = zip(d, R, y)
    cm = countmap(tuples)
    d = [key[1] for key in keys(cm)]
    R = [key[2] for key in keys(cm)]
    y = [key[3] for key in keys(cm)]
    counts = collect(values(cm))
    for index in eachindex(d, R, y, counts)
        @assert cm[(d[index], R[index], y[index])] == counts[index]
    end 
    cm2 = Dict((d[index], R[index], y[index]) => counts[index] for index in eachindex(d, R, y, counts))
    @assert cm == cm2 

    d = SizedVector{length(d), Int}(d)
    R = SizedVector{length(R), Int}(R)
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax) # contains the matrices with values y^k(1-y)^{2n-k}.
    binomialcoefficients = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n] = preparedata(n, y)
    end
    #p = Vector{Any}(undef, nmax)
    #binomialprobs = Vector{Any}(undef, nmax)
    hvals = [item for item in hvals] # change type of hvals to more specific one then Any.
    binomialcoefficients = [item for item in binomialcoefficients] # same.
    
    results = Vector{Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, BitMatrix, Matrix{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = [[MMatrix{2n+1, R+1, Float64}(undef) for R in 1:Rmax] for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Matrix{Float64}(undef, nsteps, nmax)
        τAsample = Matrix{Float64}(undef, nsteps, nmax)
        ϵsample = Matrix{Float64}(undef, nsteps, nmax)
        acceptedvec = BitMatrix(undef, nsteps, nmax+1) # BitMatrix is more efficient than Matrix{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Matrix{Float64}(undef, nsteps, nmax+1)
        # Initialise.
        τCprevious = Vector{Float64}(undef, nmax)
        τAprevious = Vector{Float64}(undef, nmax)
        ϵprevious = Vector{Float64}(undef, nmax)
        if init === nothing # draw initial value from the prior.
            nprevious = rand(prioronn)
            for index in eachindex(τCprevious, τAprevious)
                τCprevious[index], τAprevious[index] = rand(prioronτCτA)
                ϵprevious[index] = rand(prioronϵ)
            end 
        else 
            nprevious = init[1]
            for index in eachindex(τCprevious, τAprevious)
                _, τCprevious[index], τAprevious[index], ϵprevious[index] = init
            end
        end 
        nsample[1] = nprevious
        τCsample[1,:] = τCprevious
        τAsample[1,:] = τAprevious 
        ϵsample[1,:] = ϵprevious
        acceptedvec[1,:] .= true 
        # Calculate log likelihood and log prior density with initial values.
        loglikprevious = Vector{Float64}(undef, nmax)
        loglikproposal = Vector{Float64}(undef, nmax)
        for n in 1:nmax
            loglikprevious[n] = logprobderivedreads!(p[n], q[n], d, R, counts, Rmax, τCprevious[n], τAprevious[n], ϵprevious[n], n, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
        end
        # check if the initial value is valid. 
        any(loglikprevious .== -Inf) && throw(DomainError("The data does not support the initial value."))
        logpdfnprev = logpdf(prioronn, nprevious)
        logpdfpriorprevious = Vector{Float64}(undef, nmax)
        for n in 1:nmax
            logpdfpriorprevious[n] = logpdf(prioronτCτA, [τCprevious[n], τAprevious[n]]) + logpdf(prioronϵ, ϵprevious[n])
        end
        any(logpdfpriorprevious .== -Inf) && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1,2:end] .= loglikprevious .+ logpdfpriorprevious
        logjointprob[1,1] = loglikprevious[nprevious] + logpdfpriorprevious[nprevious] + logpdf(prioronn, nprevious)
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = repeat([MvNormal(Σ)], nmax)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.001
        distproposalϵ = repeat([Normal(0.0, σ)],nmax)
        scalefactor = repeat([1.8881], nmax) # = 2.38^2/4, as we have dimension of prior equal to 4.
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
                loglikproposal[n] = logprobderivedreads!(p[n], q[n], d, R, counts, Rmax, τCproposal[n], τAproposal[n], ϵproposal[n], n, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
            end 
            for n in 1:nmax
                logpdfpriorproposal[n] = logpdf(prioronτCτA, [τCproposal[n], τAproposal[n]]) + logpdf(prioronϵ, ϵproposal[n])
            end 
            logacceptanceprobability .= loglikproposal .- loglikprevious .+ logpdfpriorproposal .- logpdfpriorprevious
            acceptanceprobability .= exp.(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            for index in eachindex(acceptanceprobability)
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
            nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            #logpdfnprev = logpdf(prioronn, nprevious) # dubbele berekening.
            if nproposal == nprevious
                nsample[step] = nproposal
                acceptedvec[step, 1] = true
                logjointprob[step, 1] = logjointprob[step, nproposal+1] + logpdfnprev
            else 
                logpdfnprop = logpdf(prioronn, nproposal)
                logacceptprobn = logjointprob[step, nproposal+1] - logjointprob[step, nprevious+1] + logpdfnprop - logpdfnprev
                acceptprobn = exp(logacceptprobn)
                if rand() ≤ acceptprobn
                    nsample[step] = nproposal
                    nprevious = nproposal
                    logpdfnprev = logpdf(prioronn, nprevious)
                    acceptedvec[step, 1] = true
                    logjointprob[step, 1] = logjointprob[step, nproposal+1] + logpdfnprop
                else
                    nsample[step] = nprevious
                    acceptedvec[step, 1] = false 
                    logjointprob[step, 1] = logjointprob[step, nprevious+1] + logpdfnprev
                end 
            end             
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
                for n in 1:nmax
                    if step < nsteps÷4 && step > 1000  
                        acceptancerate = mean(acceptedvec[step-999:step, n+1]) # acceptance rate in the last 1000 steps.
                        if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                            scalefactor = 3*scalefactor/4
                            println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                        elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                            scalefactor = 4*scalefactor/3
                            println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                        end # So do nothing when the acceptancerate is between 0.1 and 0.6.     
                    end 
                    # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                    Σ = scalefactor[n] * (cov([τCsample[1:10:step, n] τAsample[1:10:step, n]], corrected=true) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                    # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                    distproposalτCτA[n] = MvNormal(Σ)
                    # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                    # so not really necessary. 
                    σ = std(ϵsample[1:10:step, n], corrected=true) + 1.0e-10 # small threshold to prevent getting zero variance. 
                    distproposalϵ[n] = Normal(0.0, sqrt(scalefactor[n]) * σ)
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
