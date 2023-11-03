"""
    module PosteriorSampling

This package provides functions to sample from the posterior. 
"""
module PosteriorSampling

# Idea for the prior on \eps. Given epsPrev, take uniform on [epsPrev/2, 3epsPrev/2]. But then it is not symmetric anymore, so q(x | y) \neq q(y | x).  

using Distributions, # for probability distributions
      Likelihood, # for the likelihood
      LinearAlgebra, # Needed for the identity matrix I. 
      #PrecompileTools, # for precompilation, to obtain faster startup times.
      StaticArrays, # for SizedVector and MVector
      Base.Threads # for multithreading.

export MHsampler, # The mcmc sampler.
       allliks, # Calclute all likelihoods in a given range of parameters.
       maxlik, # The maximum likelihood estimator.
       REXMCMC, # Replica eXchange MCMC. https://www.tweag.io/blog/2020-10-28-mcmc-intro-4/
       MCMCSimProp

# fixed coverage, fixed n, no errors. 
"""
    MHsampler(nchains::Int, nsteps::Int, prioronτCτA::ContinuousMultivariateDistribution, R::Int, d::AbstractVector{Int}, y::AbstractVector{Float64}, fixedn::Int, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64}}=nothing, temperature::Float64 = 1.0)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on τC and τA, the number of individuals is fixed. Instead of using the recommended scaling factor 2.8, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress. 

It returns a vector of length nchains, where each element are the results of of an mcmc chain. Each vector element is a 4-tuple consisting of 4 vectors: a vector with τC values, a vector with τA values, a vector with telling whether a proposal was accepted or not and the unnormalised log posterior densities. 

For example, 
using Statistics 
results = MHsampler(nchains, nsteps, prioronτCτA, R, d, y, fixedn, messages)
τCsample, τAsample, acceptedvec, logjointprob = results[1]
gives the τC, τA samples from the first chain. 

This sampler assumes that there are no errors, that the number of individuals is known, and that there is constant coverage. 
"""
function MHsampler(nchains::Int, nsteps::Int, prioronτCτA::ContinuousMultivariateDistribution, R::Int, d::AbstractVector{Int}, y::AbstractVector{Float64}, fixedn::Int, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64}}=nothing, temperature::Float64 = 1.0)
    # Check validity of the data.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads are non-negative."))
    maximum(d) ≤ R || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0"))
    length(d) == length(y) || throw(DimensionMismatch("d and y should have the same length."))
    fixedn ≥ 1 || throw(DomainError("n should be a positive integer. There is always at least one individual."))
    fixedn ≤ 10 || throw(ArgumentError("The current implementation can handle at most 10 individuals."))
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))
    # Filter out frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
    # Allocating vectors for the likelihood function.
    Q, Qꜜ, hvals, binomialcoefficients, binomialprobs = preparedata(fixedn, R, y)
    # Allocate vectors in which the MCMC samples will be stored.
    results = Vector{Tuple{Vector{Float64}, Vector{Float64}, BitVector, Vector{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results)
        p = MVector{2fixedn+1, Float64}(undef)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps)
        logjointprob = Vector{Float64}(undef, nsteps)
        # initialise the first sample.
        # First value is drawn from the prior.
        if init === nothing
            τCprevious, τAprevious = rand(prioronτCτA)
        else 
            τCprevious, τAprevious = init 
        end 
        τCsample[1] = τCprevious 
        τAsample[1] = τAprevious 
        acceptedvec[1] = true 
        # calculate log likelihood with initial values. 
        loglikprevious = logprobderivedreads!(p, d, τCprevious, τAprevious, Q, Qꜜ, binomialcoefficients, binomialprobs, hvals)
        # Calculate the log density with initial values. 
        loglikprevious == -Inf && throw(DomainError("The initial value is outside the support of the likelihood."))
        logpdfpriorprevious = logpdf(prioronτCτA, [τCprevious, τAprevious])
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
        # Initialise the proposal density. 
        scalefactor = 2.8 # ≈ 2.38^2/2, as the prior has dimension 2.
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposal = MvNormal(Σ) # The multivariate normal distribution with mean zero and covariance matrix Σ. 
        for step in 2:nsteps 
            # sample proposals.
            τCproposal, τAproposal = [τCprevious, τAprevious] + rand(distproposal)
            loglikproposal = logprobderivedreads!(p, d, τCproposal, τAproposal, Q, Qꜜ, binomialcoefficients, binomialprobs, hvals)
            logpdfpriorproposal = logpdf(prioronτCτA, [τCproposal, τAproposal])
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + logpdfpriorproposal - logpdfpriorprevious )
            acceptanceprobability = exp(logacceptanceprobability)
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal ) 
                # update the values for the next iteration.
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
                # loglikprevious and logpdfpriorprevious remain unchanged.
            end 
            if step % 1000 == 0
                if step < nsteps÷4 && step > 1000
                    acceptancerate = mean(acceptedvec[step-999:step])
                    if acceptancerate < 0.1
                        scalefactor = 3*scalefactor/4
                        println("I, chain ", chainindex, " talking, ", " I scaled scalefactor with 3/4, new scalefactor = ", scalefactor, '.')
                    elseif acceptancerate > 0.6
                        scalefactor = 4*scalefactor/3
                        println("I, chain ", chainindex, " talking, ", " I scaled scalefactor with 4/3, new scalefactor = ", scalefactor, '.')
                    end # So do nothing when the acceptancerate is between 0.2 and 0.5. 
                end # if step < nsteps÷4
                # We use Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2),
                # the covariance is calculated from the previous steps and used in the next iterations.
                Σ = scalefactor * (cov([τCsample[begin:10:step] τAsample[begin:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. 
                distproposal = MvNormal(Σ)
            end
            if messages > 0 # no messages for messages ≤ 0. 
                if step % messages == 0
                    println("Chain ", chainindex, " talking, I finished step ", step, " of ", nsteps, ", ", round(Int, 100*step / nsteps), " percent finished.")
                end
            end 
        end
        results[chainindex] = (τCsample, τAsample, acceptedvec, logjointprob)
    end
    return results
end

# Fixed R, no errors. 
"""
    MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on n, τC, and τA. Instead of using the recommended scaling factor 1.888, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress.

The initial values are choosen from the prior.

This function returns a vector of length nchains. Each element contains the results of one chain; it is a tuple of five vectors: a vector containing the n-values, a vector containing the τC-values, a vector containing the τA-values, a vector telling whether a certain proposal was accepted or not and a vector with unnormalised posterior log-densities. So for example 

```julia
results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, R, d, y, messages)
nsample, τCsample, τAsample, acceptedvec, logjointprob = results[1]
```

gives the results of the first chain.
"""
function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    # De conditie hieronder kan relatief eenvoudig gefixt worden door de waarden van te voren uit te rekenen tot en met r, zodat 
    # P(X>r)^nsteps ≤ 0.1. Als we dan het punt raken waar n>r, dan moeten we de vector uitbreiden. 
    maximum(prioronn) < Inf || throw(ArgumentError("For the current implementation we require that the prior on n is bounded. I.e. there is an N, so that P(prioronn ≤ N)=1"))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    R > 0 || throw(DomainError("All coverages should be positive."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
    nmax = maximum(prioronn)
    nmax ≤ 10 || throw(ArgumentError("The current implementation can handle at most 10 individuals."))
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax)
    binomialcoefficients = Vector{Any}(undef, nmax)
    #p = Vector{Any}(undef, nmax)
    binomialprobs = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n], binomialprobs[n] = preparedata(n, R, y)
    end 
    # Q, Qꜜ, hvals, binomialcoefficients, binomialprobs = preparedata(fixedn, R, y)
    hvals = [item for item in hvals]
    binomialcoefficients = [item for item in binomialcoefficients]
    #p = [item for item in p]
    binomialprobs = [item for item in binomialprobs]
    results = Vector{Tuple{Vector{Int}, Vector{Float64}, Vector{Float64}, BitVector, Vector{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps) # BitVector is more efficient than Vector{Bool} (roughly, one bit per item, instead of one byte.)
        logjointprob = Vector{Float64}(undef, nsteps)
        # Initialise.
        if init === nothing
            nprevious = rand(prioronn)
            τCprevious, τAprevious = rand(prioronτCτA)
        else 
            nprevious, τCprevious, τAprevious = init 
        end 
        nsample[1] = nprevious
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious 
        acceptedvec[1] = true 
        # Calculate log likelihood and log prior density with initial values.
        loglikprevious = logprobderivedreads!(p[nprevious], d, τCprevious, τAprevious, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], binomialprobs[nprevious], hvals[nprevious])
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronn, nprevious) + logpdf(prioronτCτA, [τCprevious, τAprevious])
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = MvNormal(Σ)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        scalefactor = 1.888 # 1.888 ≈ 2.38^2/3, as we have dimension of prior equal to 3.
        #distproposaln = Categorical([1/3, 1/3, 1/3]) # takes values 1, 2 and 3 with probabilities 0.1, 0.8 and 0.1. 
        for step in 2:nsteps  
            # draw proposals. 
            τCproposal, τAproposal = [τCprevious, τAprevious] + rand(distproposalτCτA)
            # The distribution below is symmetric, i.e. P(a;b)=P(b;a), even for a,b equal or near 1 or n.  
            nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            nproposal = rand(Categorical(nproposalprobvec))
            # calculate the loglikelihood
            loglikproposal = logprobderivedreads!(p[nproposal], d, τCproposal, τAproposal, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], binomialprobs[nproposal], hvals[nproposal])
            logpdfpriorproposal = logpdf(prioronn, nproposal) + logpdf(prioronτCτA, [τCproposal, τAproposal])
            # calculate the acceptance probability. 
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + logpdfpriorproposal - logpdfpriorprevious )
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution on the interval [0,1].
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                nsample[step] = nprevious = nproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal )
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                nsample[step] = nprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
                # no updating: loglikprevious and logpdfpriorprevious remain the same.
            end 
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
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
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
            end
            if messages > 0 # No messages at all when messages ≤ 0.
                if step % messages == 0
                    println("Chain ", chainindex, " talking, I finished step ", step, " of ", nsteps, ", ", round(Int, 100*step / nsteps, RoundDown), " percent finished.")
                end
            end
        end # MCMC step
        results[chainindex] = (nsample, τCsample, τAsample, acceptedvec, logjointprob)
    end 
    return results
end

# fixed R, error model, unknown n.
# init should also provide a value for ϵ. DONE
# to do: instantiate q. DONE
"""
    MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)

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
function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)
    # Check conditions.
    length(prioronτCτA) == 2 || throw(DimensionMismatch("prior on [τC, τA] must be of length 2."))
    minimum(prioronn) ≥ 1 || throw(ArgumentError("prioronn should have support on positive integers. There is always at least one individual."))
    maximum(prioronn) < Inf || throw(ArgumentError("For the current implementation we require that the prior on n is bounded. I.e. there is an N, so that P(prioronn ≤ N)=1"))
    minimum(prioronϵ) ≥ 0 || throw(ArgumentError("prioronϵ should have support on [0,1], but the minimum is negative."))
    maximum(prioronϵ) ≤ 1 || throw(ArgumentError("prioronϵ should have support on [0,1], but the maximum exceeds 1."))
    nchains > 0 || throw(ArgumentError("There should be at least one chain."))
    nsteps > 0 || throw(ArgumentError("nsteps should be positive."))
    R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    maximum(d) ≤  R || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
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
        loglikprevious = logprobderivedreads!(p[nprevious], q[nprevious], d, τCprevious, τAprevious, ϵprevious, nprevious, R, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronn, nprevious) + 
            logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
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
            loglikproposal = logprobderivedreads!(p[nproposal], q[nproposal], d, τCproposal, τAproposal, ϵproposal, nproposal, R, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])
            logpdfpriorproposal = logpdf(prioronn, nproposal) + logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious ) 
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                nsample[step] = nprevious = nproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal )
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                nsample[step] = nprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
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
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step]) + 1.0e-5 # small threshold to prevent getting zero variance. 
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

# variable R, error model. 
# logprobderivedreads!(p::AbstractVector{Float64}, q::AbstractVector{<:AbstractMatrix{Float64}}, reads::AbstractVector{Int}, coverages::AbstractVector{Int}, maxcoverage::Int, τC::Float64, τA::Float64, ϵ::Float64, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
"""
    MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on n, τC, τA and ϵ. Instead of using the recommended scaling factor 1.4161, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress.

If not specified, the initial values are choosen from the prior.

This function returns a vector of length nchains. Each element contains the results of one chain; it is a tuple of six vectors: a vector containing the n-values, a vector containing the τC-values, a vector containing the τA-values, a vector containing the ϵ values, a vector telling whether a certain proposal was accepted or not and a vector with unnormalised posterior log-densities. So for example 

```julia
results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, R, d, y, messages)
nsample, τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[1]
```

gives the results of the first chain.
"""
function MHsampler(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)
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
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    R = R[allowedindices]
    R = SizedVector{length(R), Int}(R)
    y = y[allowedindices]
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
        loglikprevious = logprobderivedreads!(p[nprevious], q[nprevious], d, R, Rmax, τCprevious, τAprevious, ϵprevious, nprevious, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronn, nprevious) + 
            logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
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
            loglikproposal = logprobderivedreads!(p[nproposal], q[nproposal], d, R, Rmax, τCproposal, τAproposal, ϵproposal, nproposal, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])
            logpdfpriorproposal = logpdf(prioronn, nproposal) + logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious ) 
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                nsample[step] = nprevious = nproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal )
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                nsample[step] = nprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
                # no updating: loglikprevious and logpdfpriorprevious remain the same.
            end 
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
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
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step]) + 1.0e-5 # small threshold to prevent getting zero variance. 
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
## variable R, error model, fixed n. 
##
##
"""
    MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)

Adaptive Metropolis (Haario et al. 2001, mentioned in Rosenthal 2010 "Optimal proposal distributions and adaptive MCMC, section 3.2), with prior on τC, τA and ϵ, the true number of individuals is assumed to be known. Instead of using the recommended scaling factor 1.4161, we use the first nsteps÷4 steps to find a scaling factor that results in an acceptance ratio between 0.1 and 0.6. Setting messages ≤ 0, results in no messages being print, if message>0, then every message steps, a message is print, updating the progress.

If not specified, the initial values are choosen from the prior.

This function returns a vector of length nchains. Each element contains the results of one chain; it is a tuple of five vectors: a vector containing the a vector containing the τC-values, a vector containing the τA-values, a vector containing the ϵ values, a vector telling whether a certain proposal was accepted or not and a vector with unnormalised posterior log-densities. So for example 

```julia
results = MHsampler(nchains, nsteps, prioronn, prioronτCτA, R, d, y, messages)
τCsample, τAsample, ϵsample, acceptedvec, logjointprob = results[1]
```

gives the results of the first chain.

Here we assume that the true number of individuals in the sediment is known.
"""
function MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)
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
    Rmax = maximum(R)
    Rmax > 0 || throw(ArgumentError("Inference is impossible when all coverages are zero."))
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    R = R[allowedindices]
    R = SizedVector{length(R), Int}(R)
    y = y[allowedindices]
    Q, Qꜜ, hvals, binomialcoefficients = preparedata(fixedn, y)
    results = Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, BitVector, Vector{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = MVector{2fixedn+1, Float64}(undef)
        q = [MMatrix{2fixedn+1, R+1, Float64}(undef) for R in 1:Rmax]
        #nsample = Vector{Int}(undef, nsteps)
        τCsample = Vector{Float64}(undef, nsteps)
        τAsample = Vector{Float64}(undef, nsteps)
        ϵsample = Vector{Float64}(undef, nsteps)
        acceptedvec = BitVector(undef, nsteps) # BitVector is more efficient than Vector{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Vector{Float64}(undef, nsteps)
        # Initialise.
        if init === nothing # draw initial value from the prior.
            τCprevious, τAprevious = rand(prioronτCτA)
            ϵprevious = rand(prioronϵ)
        else
            #nprevious, 
            τCprevious, τAprevious, ϵprevious = init 
        end
        τCsample[1] = τCprevious
        τAsample[1] = τAprevious
        ϵsample[1] = ϵprevious
        acceptedvec[1] = true
        loglikprevious = logprobderivedreads!(p, q, d, R, Rmax, τCprevious, τAprevious, ϵprevious, fixedn, Q, Qꜜ, binomialcoefficients, hvals)
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
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
            #nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            #nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            ϵproposal = ϵprevious + rand(distproposalϵ)
            # calculate the loglikelihood
            loglikproposal = logprobderivedreads!(p, q, d, R, Rmax, τCproposal, τAproposal, ϵproposal, fixedn, Q, Qꜜ, binomialcoefficients, hvals)
            logpdfpriorproposal = logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious ) 
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                #nsample[step] = nprevious = nproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal )
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                #nsample[step] = nprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
                # no updating: loglikprevious and logpdfpriorprevious remain the same.
            end 
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
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
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step]) + 1.0e-5 # small threshold to prevent getting zero variance. 
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
# End variable R, error model. 
##
##
## End variable R, error model, fixed n. 
##
##

##
##
## fixed R, error model, fixed n. 
##
##
"""
    MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)

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
function MHsampler(nchains::Int, nsteps::Int, fixedn::Int, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)
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
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    #R = R[allowedindices]
    #R = SizedVector{length(R), Int}(R)
    y = y[allowedindices]
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
        loglikprevious = logprobderivedreads!(p, q, d, τCprevious, τAprevious, ϵprevious, fixedn, R, Q, Qꜜ, binomialcoefficients, hvals)
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = MvNormal(Σ)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.001
        distproposalϵ = Normal(0.0, σ)
        scalefactor = 1.8881 # = 2.38^2/3, as we have dimension of prior equal to 3.
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
            loglikproposal = logprobderivedreads!(p, q, d, τCproposal, τAproposal, ϵproposal, fixedn, R, Q, Qꜜ, binomialcoefficients, hvals)
            logpdfpriorproposal = logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious ) 
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal )
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
                # no updating: loglikprevious and logpdfpriorprevious remain the same.
            end 
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
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
                # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step]) + 1.0e-5 # small threshold to prevent getting zero variance. 
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
##


##
##
## MCMC with simulated proposals. 
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
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
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
    results = Vector{Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, BitMatrix, Matrix{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Matrix{Float64}(undef, nsteps, nmax)
        τAsample = Matrix{Float64}(undef, nsteps, nmax)
        ϵsample = Matrix{Float64}(undef, nsteps, nmax)
        acceptedvec = BitMatrix(undef, nsteps, nmax+1) # BitVector is more efficient than Vector{Bool} 
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
            loglikprevious[n] = logprobderivedreads!(p[n], q[n], d, τCprevious[n], τAprevious[n], ϵprevious[n], n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
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
                loglikproposal[n] = logprobderivedreads!(p[n], q[n], d, τCproposal[n], τAproposal[n], ϵproposal[n], n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
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
                logjointprob[step, 1] = logjointprob[step, nprevious] + logpdfnprev
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
                            scalefactor[n] = 3*scalefactor[n]/4
                            #println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                        elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                            scalefactor[n] = 4*scalefactor[n]/3
                            #println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                        end # So do nothing when the acceptancerate is between 0.1 and 0.15. 
                    end 
                    # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                    Σ = scalefactor[n] * (cov([τCsample[1:10:step, n] τAsample[1:10:step, n]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                    # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                    distproposalτCτA[n] = MvNormal(Σ)
                    # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                    # so not really necessary. 
                    σ = std(ϵsample[1:10:step, n]) + 1.0e-10 # small threshold to prevent getting zero variance. 
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
## end new, adding error model.
##
##
##


#= backup 
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
    
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
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
    results = Vector{Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, BitMatrix, Matrix{Float64}}}(undef,nchains)
    @threads for chainindex in eachindex(results) 
        p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
        q = [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
        nsample = Vector{Int}(undef, nsteps)
        τCsample = Matrix{Float64}(undef, nsteps, nmax+1)
        τAsample = Matrix{Float64}(undef, nsteps, nmax+1)
        ϵsample = Matrix{Float64}(undef, nsteps, nmax+1)
        acceptedvec = BitMatrix(undef, nsteps, nmax+1) # BitVector is more efficient than Vector{Bool} 
        # (roughly, one bit per item, instead of one byte.)
        logjointprob = Matrix{Float64}(undef, nsteps, nmax+1)
        # Initialise.
        τCprevious = Vector{Float64}(undef, nmax+1)
        τAprevious = Vector{Float64}(undef, nmax+1)
        ϵprevious = Vector{Float64}(undef, nmax+1)
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
        loglikprevious = Vector{Float64}(undef, nmax+1)
        loglikproposal = Vector{Float64}(undef, nmax+1)
        loglikprevious[1] = logprobderivedreads!(p[nprevious], q[nprevious], d, τCprevious[1], τAprevious[1], ϵprevious[1], nprevious, R, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])
        for n in 1:nmax
            loglikprevious[n+1] = logprobderivedreads!(p[n], q[n], d, τCprevious[n+1], τAprevious[n+1], ϵprevious[n+1], n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
        end
        # check if the initial value is valid. 
        any(loglikprevious .== -Inf) && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = Vector{Float64}(undef, nmax+1)
        logpdfpriorprevious[1] = logpdf(prioronn, nprevious) + 
            logpdf(prioronτCτA, [τCprevious[1], τAprevious[1]]) + logpdf(prioronϵ, ϵprevious[1])
        for n in 1:nmax
            logpdfpriorprevious[n+1] = logpdf(prioronτCτA, [τCprevious[n+1], τAprevious[n+1]]) + logpdf(prioronϵ, ϵprevious[n+1])
        end 
        any(logpdfpriorprevious .== -Inf) && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1,:] = loglikprevious + logpdfpriorprevious
        # Define proposal densities. 
        Σ = [0.00001 0.0; 0.0 0.00001]
        distproposalτCτA = repeat([MvNormal(Σ)], nmax+1)# Multivariate normal distribution with mean zero and 
        # covariance matrix Σ.
        σ = 0.001
        distproposalϵ = repeat([Normal(0.0, σ)],nmax+1)
        scalefactor = vcat([1.4161], repeat([1.8881], nmax)) # = 2.38^2/4, as we have dimension of prior equal to 4.
        τCproposal = Vector{Float64}(undef, nmax+1)
        τAproposal = Vector{Float64}(undef, nmax+1)
        ϵproposal = Vector{Float64}(undef, nmax+1)
        logpdfpriorproposal = Vector{Float64}(undef, nmax+1)
        logacceptanceprobability = Vector{Float64}(undef, nmax+1)
        acceptanceprobability = Vector{Float64}(undef, nmax+1)
        accepted = BitVector(undef, nmax+1)
        for step in 2:nsteps  
            # draw proposals.
            nproposalprobvec = vcat(repeat([0.2 / (nmax-1)], nprevious-1), [0.8], repeat([0.2 / (nmax-1)], nmax - nprevious))
            nproposal = rand(Categorical(nproposalprobvec)) # so with 80% chance nproposal is 
            # equal to nprevious. 
            if nproposal == nprevious
                for index in eachindex(τCproposal, τAproposal) 
                    τCproposal[index], τAproposal[index] = [τCprevious[index], τAprevious[index]] + rand(distproposalτCτA[index])
                end
            else
                τCproposal[1], τAproposal[1] = [τCprevious[nproposal+1], τAprevious[nproposal+1]] + rand(distproposalτCτA[1])
                for index in 2:nmax+1 
                    τCproposal[index], τAproposal[index] = [τCprevious[index], τAprevious[index]] + rand(distproposalτCτA[index])
                end
            end 
            for index in eachindex(ϵproposal, distproposalϵ) 
                ϵproposal[index] = ϵprevious[index] + rand(distproposalϵ[index])
            end
            # calculate the loglikelihood
            loglikproposal[1] = logprobderivedreads!(p[nproposal], q[nproposal], d, τCproposal[1], τAproposal[1], ϵproposal[1], nproposal, R, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])
            for n in 1:nmax
                loglikproposal[n+1] = logprobderivedreads!(p[n], q[n], d, τCproposal[n+1], τAproposal[n+1], ϵproposal[n+1], n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
            end 
            logpdfpriorproposal[1] = logpdf(prioronn, nproposal) + logpdf(prioronτCτA, [τCproposal[1], τAproposal[1]]) + logpdf(prioronϵ, ϵproposal[1])
            for n in 1:nmax
                logpdfpriorproposal[n+1] = logpdf(prioronτCτA, [τCproposal[n+1], τAproposal[n+1]]) + logpdf(prioronϵ, ϵproposal[n+1])
            end 
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability .= loglikproposal .- loglikprevious .+ logpdfpriorproposal .- logpdfpriorprevious
            acceptanceprobability .= exp.(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            for index in eachindex(acceptanceprobability)
                accepted[index] = rand() ≤ acceptanceprobability[index]
            end 
            if accepted[1]
                nsample[step] = nprevious = nproposal
            else
                nsample[step] = nprevious
            end
            for index in eachindex(accepted)
                if accepted[index]
                    τCsample[step, index] = τCprevious[index] = τCproposal[index]
                    τAsample[step, index] = τAprevious[index] = τAproposal[index]
                    ϵsample[step, index] = ϵprevious[index] = ϵproposal[index]
                    acceptedvec[step, index] = true
                    logjointprob[step, index] = loglikproposal[index] + logpdfpriorproposal[index]
                    # Update everything for the next cycle. 
                    loglikprevious[index] = loglikproposal[index]
                    logpdfpriorprevious[index] = logpdfpriorproposal[index]
                else
                    τCsample[step, index] = τCprevious[index]
                    τAsample[step, index] = τAprevious[index]
                    ϵsample[step, index] = ϵprevious[index]
                    acceptedvec[step, index] = false
                    logjointprob[step, index] = loglikprevious[index] + logpdfpriorprevious[index]
                    # no updating: loglikprevious and logpdfpriorprevious remain the same.
                end 
            end 
            if step % 1000 == 0
                # the covariance is calculated from the previous steps and used in the next iteration.
                for index in eachindex(scalefactor)
                    if step < nsteps÷4 && step > 1000  
                        acceptancerate = mean(acceptedvec[step-999:step, index]) # acceptance rate in the last 1000 steps.
                        if acceptancerate < 0.1 # acceptancerate too low, take smaller steps.
                            scalefactor[index] = 3*scalefactor[index]/4
                            #println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 3/4, new scalefactor = ", scalefactor, '.')
                        elseif acceptancerate > 0.15 # acceptancerate too large, take bigger steps.
                            scalefactor[index] = 4*scalefactor[index]/3
                            #println("Chain ", chainindex, " talking, I NEW scaled scalefactor with factor 4/3, new scalefactor = ", scalefactor, '.')
                        end # So do nothing when the acceptancerate is between 0.1 and 0.6. 
                    end 
                    # update Σ to the empirical covariance matrix of [τC, τA] based on the first step steps. 
                    Σ = scalefactor[index] * (cov([τCsample[1:10:step, index] τAsample[1:10:step, index]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                    # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                    distproposalτCτA[index] = MvNormal(Σ)
                    # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                    # so not really necessary. 
                    σ = std(ϵsample[1:10:step, index]) + 1.0e-5 # small threshold to prevent getting zero variance. 
                    distproposalϵ[index] = Normal(0.0, sqrt(scalefactor[index]) * σ)
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
=#


##
##
## Implementation of replica exchange MCMC
##
##
"""
    REXMCMC(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing, temperature::Float64 = 1.0)

TBW
"""
function REXMCMC(nchains::Int, nsteps::Int, prioronn::DiscreteUnivariateDistribution, prioronτCτA::ContinuousMultivariateDistribution, prioronϵ::ContinuousUnivariateDistribution, R::Int, d::Vector{Int}, y::Vector{Float64}, αs::Vector{Float64}, messages::Int=nsteps÷100; init::Union{Nothing, Tuple{Int64, Float64, Float64, Float64}}=nothing)
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
    0.0 < temperature ≤ 1.0 || throw(DomainError("The temperature should satisfy 0.0 < temperature ≤ 1.0."))
    αs[1] == 1.0 || throw(DomainError("αs[1] should be 1.0"))
    all(0.0 .< αs .≤ 1.0) || throw(DomainError("All αs should be strictly larger than 0.0 and smaller or equal to 1.0."))
    for i in 1:length(αs)-1
        αs[i] > αs[i+1] || throw(DomainError("αs[$i] is not strictly larger than αs[$(i+1)]."))
    end 
    ##
    ##
    ## HIER VERDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ##
    ##
    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
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
        loglikprevious = logprobderivedreads!(p[nprevious], q[nprevious], d, τCprevious, τAprevious, ϵprevious, nprevious, R, Qs[nprevious], Qꜜs[nprevious], binomialcoefficients[nprevious], hvals[nprevious])
        # check if the initial value is valid. 
        loglikprevious == -Inf && throw(DomainError("The data does not support the initial value."))
        logpdfpriorprevious = logpdf(prioronn, nprevious) + 
            logpdf(prioronτCτA, [τCprevious, τAprevious]) + logpdf(prioronϵ, ϵprevious)
        logpdfpriorprevious == -Inf && throw(DomainError("The initial value is outside the support of the prior."))
        logjointprob[1] = temperature * ( loglikprevious + logpdfpriorprevious )
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
            loglikproposal = logprobderivedreads!(p[nproposal], q[nproposal], d, τCproposal, τAproposal, ϵproposal, nproposal, R, Qs[nproposal], Qꜜs[nproposal], binomialcoefficients[nproposal], hvals[nproposal])
            logpdfpriorproposal = logpdf(prioronn, nproposal) + logpdf(prioronτCτA, [τCproposal, τAproposal]) + logpdf(prioronϵ, ϵproposal)
            # calculate the acceptance probability. 
            # Note that our proposal function for n, τC, τA, ϵ is symmetric, i.e. q(x|y)=q(y|x).
            logacceptanceprobability = temperature * ( loglikproposal - loglikprevious + 
                logpdfpriorproposal - logpdfpriorprevious ) 
            acceptanceprobability = exp(logacceptanceprobability)
            # accept or reject the proposal with probability max(1, acceptanceprobability), note that rand() samples from the uniform distribution.
            accepted = rand() ≤ acceptanceprobability
            if accepted
                τCsample[step] = τCprevious = τCproposal
                τAsample[step] = τAprevious = τAproposal
                nsample[step] = nprevious = nproposal
                ϵsample[step] = ϵprevious = ϵproposal
                acceptedvec[step] = true
                logjointprob[step] = temperature * ( loglikproposal + logpdfpriorproposal )
                # Update everything for the next cycle. 
                loglikprevious = loglikproposal
                logpdfpriorprevious = logpdfpriorproposal
            else
                τCsample[step] = τCprevious
                τAsample[step] = τAprevious
                nsample[step] = nprevious
                ϵsample[step] = ϵprevious
                acceptedvec[step] = false
                logjointprob[step] = temperature * ( loglikprevious + logpdfpriorprevious )
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
                Σ = scalefactor * (cov([τCsample[1:10:step] τAsample[1:10:step]]) + 1.0e-10*I) # 1.0e-10*I is a small treshold, to prevent getting (almost) the zero matrix. In addition it guarantees that Σ is positive definite, as the empirical covariance matrix is positive 
                # semi-definite, so adding a positive definite matrix makes the whole positive definite. 
                distproposalτCτA = MvNormal(Σ)
                # Do I need to do something for the discrete parameter? NO, in the simulations, the sampler gets stuck in one value of n,
                # so not really necessary. 
                σ = std(ϵsample[1:10:step]) + 1.0e-5 # small threshold to prevent getting zero variance. 
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
##
##
## End implementation of replica exchange MCMC. 
##
##

# without errors. Fixed R. 
"""
    allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, R::Int, d::Vector{Int}, y::Vector{Float64})

Calclute all log likelihoods in the ranges of n, τC and τA values. Return a dictionary with keys tuples (n, τC, τA) and as value the corresponding log likelihood.

Here we assume that there is no error and that the coverage is fixed.
"""
function allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, R::Int, d::Vector{Int}, y::Vector{Float64})
    # Check conditions.
    R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
    nmax = maximum(nrange)
    # Misschien kan ik het symmetrie probleem oplossen door nmax+1 te nemen. Dan is de prior probibility 0 en krijg je sowieso rejection.
    # maar het probleem is dan nog steeds niet opgelost voor n = 0. Je moet ook kijken wat het numeriek betekend.
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax)
    binomialcoefficients = Vector{Any}(undef, nmax)
    p = Vector{Any}(undef, nmax)
    binomialprobs = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n], binomialprobs[n] = preparedata(n, R, y)
        p[n] = MVector{2n+1, Float64}(undef)
    end 
    hvals = [item for item in hvals]
    binomialcoefficients = [item for item in binomialcoefficients]
    p = [item for item in p]
    binomialprobs = [item for item in binomialprobs]
    logliks = Dict{Tuple{Int64, Float64, Float64}, Float64}()
    for n in nrange
        for τC in τCrange
            for τA in τArange
                logliks[(n, τC, τA)] = logprobderivedreads!(p[n], d, τC, τA, Qs[n], Qꜜs[n], binomialcoefficients[n], binomialprobs[n], hvals[n])
            end
        end 
    end 
    return logliks
end

# with errors. Fixed R. 
"""
    allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, ϵrange::AbstractVector{Float64}, R::Int, d::Vector{Int}, y::Vector{Float64})

Calclute all log likelihoods in the ranges of n, τC, τA and ϵ values. Return a dictionary with keys tuples (n, τC, τA, ϵ) and as value the corresponding log likelihood.

Here we assume that the coverage is fixed
"""
function allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, ϵrange::AbstractVector{Float64}, R::Int, d::Vector{Int}, y::Vector{Float64})
    # Check conditions.
    R > 0 || throw(DomainError("All coverages should be positive.")) # Are zero coverages allowed?
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))

    # Remove frequencies 0.0 and 1.0
    allowedindices = 0.0 .< y .< 1.0
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    y = y[allowedindices]
    nmax = maximum(nrange)
    # Misschien kan ik het symmetrie probleem oplossen door nmax+1 te nemen. Dan is de prior probibility 0 en krijg je sowieso rejection.
    # maar het probleem is dan nog steeds niet opgelost voor n = 0. Je moet ook kijken wat het numeriek betekend.

    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax) # contains the matrices with values y^k(1-y)^{2n-k}.
    binomialcoefficients = Vector{Any}(undef, nmax)
    #p = Vector{Any}(undef, nmax)
    #binomialprobs = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n], _ = preparedata(n, R, y)
    end
    hvals = [item for item in hvals] # change type of hvals to more specific one then Any.
    binomialcoefficients = [item for item in binomialcoefficients] # same.
    p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
    q = [MMatrix{2n+1, R+1, Float64}(undef) for n in 1:nmax]
    logliks = Dict{Tuple{Int64, Float64, Float64, Float64}, Float64}()
    for n in nrange
        for τC in τCrange
            for τA in τArange
                for ϵ in ϵrange
                    logliks[(n, τC, τA, ϵ)] = logprobderivedreads!(p[n], q[n], d, τC, τA, ϵ, n, R, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
                end  
            end 
        end 
    end 
    return logliks 
end

# with errors. Flexible R. 
"""
    allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, ϵrange::AbstractVector{Float64}, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64})

Calclute all log likelihoods in the ranges of n, τC, τA and ϵ values. Return a dictionary with keys tuples (n, τC, τA, ϵ) and as value the corresponding log likelihood.
"""
function allliks(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, ϵrange::AbstractVector{Float64}, R::Vector{Int}, d::Vector{Int}, y::Vector{Float64})
    # Check conditions.
    minimum(R) ≥ 0 || throw(DomainError("All coverages should be non-negative.")) # Are zero coverages allowed?
    Rmax = maximum(R)
    Rmax > 0 || throw(DomainError("Inference is impossible with all coverages zero."))
    minimum(d) ≥ 0 || throw(DomainError("All reads should be non-negative."))
    all(d .≤ R) || throw(DomainError("All reads should be at most the coverage."))
    all(0.0 .≤ y .≤ 1.0) || throw(DomainError("The frequencies should be between 0.0 and 1.0."))
    length(d) == length(y) || throw(DimensionMismatch("d and y should be of the same length."))
    # Remove frequencies 0.0 and 1.0
    allowedindices = (0.0 .< y .< 1.0) .& (R .> 0)
    d = d[allowedindices]
    d = SizedVector{length(d), Int}(d)
    R = R[allowedindices]
    R = SizedVector{length(R), Int}(R)
    y = y[allowedindices]
    nmax = maximum(nrange)
    Qs = Vector{EigenExpansion}(undef, nmax)
    Qꜜs = Vector{EigenExpansion}(undef, nmax)
    hvals = Vector{Any}(undef, nmax) # contains the matrices with values y^k(1-y)^{2n-k}.
    binomialcoefficients = Vector{Any}(undef, nmax)
    for n in 1:nmax
        Qs[n], Qꜜs[n], hvals[n], binomialcoefficients[n] = preparedata(n, y)
    end
    hvals = [item for item in hvals] # change type of hvals to more specific one then Any.
    binomialcoefficients = [item for item in binomialcoefficients] # same.
    p = [MVector{2n+1, Float64}(undef) for n in 1:nmax]
    q = [[MMatrix{2n+1, R+1, Float64}(undef) for R in 1:Rmax] for n in 1:nmax]
    logliks = Dict{Tuple{Int64, Float64, Float64, Float64}, Float64}()
    for n in nrange
        for τC in τCrange
            for τA in τArange
                for ϵ in ϵrange
                    logliks[(n, τC, τA, ϵ)] = logprobderivedreads!(p[n], q[n], d, R, Rmax, τC, τA, ϵ, n, Qs[n], Qꜜs[n], binomialcoefficients[n], hvals[n])
                end  
            end 
        end 
    end 
    return logliks 
end

"""
    maxlik(args...)

Calculate the maximum likelihood estimator, with arguments from a method of allliks. 
"""
maxlik(args...) = argmax(allliks(args...))
#maxlik(nrange::AbstractVector{Int}, τCrange::AbstractVector{Float64}, τArange::AbstractVector{Float64}, R::Int, d::Vector{Int}, y::Vector{Float64}) = argmax(allliks(nrange, τCrange, τArange, R, d, y))

# Precompilation, in order to make startup time faster.
#=
@setup_workload begin
    nchains = 2
    nsteps = 3
    #init = (τC = 0.3, τA = 0.3)
    prior = Product([Uniform(), Uniform()])
    coverage = 3
    derivedreadsvec = [1,0,2]
    freqs = [0.1,0.2,0.1]
    truenumind = 1

    @compile_workload begin
        MHsampler(nchains, nsteps, prior, coverage, derivedreadsvec, freqs, truenumind, 25_000)
        MHsampler(nchains, nsteps, prior, coverage, derivedreadsvec, freqs, truenumind)  
    end
end

@setup_workload begin
    nsteps = 3
    prioronτCτA = Product([Uniform(), Uniform()])
    nchains = 2
    prioronn = DiscreteUniform(1,3)
    coverage = 3
    derivedreadsvec = [1,0,2]
    freqs = [0.1,0.2,0.1]
    @compile_workload begin
        MHsampler(nchains, nsteps, prioronn, prioronτCτA, coverage, derivedreadsvec, freqs, 25_000)
        MHsampler(nchains, nsteps, prioronn, prioronτCτA, coverage, derivedreadsvec, freqs)  
    end
end
=# 

end # package PosteriorSampling
