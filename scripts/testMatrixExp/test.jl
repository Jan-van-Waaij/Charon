using LinearAlgebra
using Likelihood: makeQ, makeQꜜ, calcmatrix
using Likelihood
using LogExpFunctions
using Random 
using StatsBase 
using Distributions 
using SpecialFunctions
using Dates 

function logprobderivedreadstestlogspace!(logq::AbstractMatrix{Float64}, intermediateresult::AbstractVector{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, logbinomialcoefficients::AbstractVector{Float64}, hvals::AbstractMatrix{Float64})
    if τC < 0 || τA < 0 || ϵ < 0 || ϵ > 1
        return -Inf 
    else
        M = calcmatrix(Q, Qꜜ, τC, τA)
        # update q. q[k+1, d+1] is the probability of d derived alleles given k derived alleles (out of 2n)
        # and error probability ϵ.
        for k in 0:2n
            dist = Binomial(R, k*(1.0-ϵ)/(2n) + (2n-k)*ϵ/(2n))
            for d in 0:R
                logq[k+1, d+1] = logbinomialcoefficients[k+1] + logpdf(dist, d)
            end
        end 
        loglik = 0.0
        for index in eachindex(reads, counts)
            mul!(intermediateresult, M, view(hvals, :, index)) # caculate M*hvals[:, index] and put it in intermediateresults.
            intermediateresult .= log.(max.(0.0, intermediateresult)) # correct for possible negative values due to numerical errors, and calculate the log. This is an in-place calculation, no new vector is created. 
            intermediateresult .+= view(logq, :, reads[index]+1) 
            loglik += counts[index] * logsumexp(intermediateresult) 
        end
        return loglik
    end 
end # Constant coverage and with errors. 

function logprobderivedreadstestexpspace!(q::AbstractArray{Float64}, reads::AbstractVector{Int}, R::Int, uniqueR::Int, counts::AbstractVector{Int}, τC::Float64, τA::Float64, ϵ::Float64, n::Int, Q::EigenExpansion, Qꜜ::EigenExpansion, binomialcoefficients::AbstractVector{Int}, hvals::AbstractMatrix{Float64})
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


function testlogspace(n, ϵ, R, expsize)
    Random.seed!(1987)
    y = rand(0.1:0.1:0.9, expsize)
    logq = Matrix{Float64}(undef, 2n+1, R+1)
    intermediateresult = Vector{Float64}(undef, 2n+1)
    Random.seed!(1987)
    reads = rand(0:R, expsize)
    cm = countmap(zip(y, reads))
    y = [key[1] for key in keys(cm)]
    reads = [key[2] for key in keys(cm)]
    counts = collect(values(cm))
    Q, Qꜜ, hvals, _ = preparedata(n, y)
    logbinomialcoefficients = [logabsbinomial(2n, k)[1] for k in 0:2n]

    logpdffunction = logprobderivedreadstestlogspace!(logq, intermediateresult, reads, R, R, counts, 0.0, 0.0, ϵ, n, Q, Qꜜ, logbinomialcoefficients, hvals)
    
    logpdftheoretical = 0.0
    logprobs = Vector{Float64}(undef, 2n+1)
    for index in eachindex(reads, counts, y)
        d = reads[index]
        distk = Binomial(2n, y[index])
        for k in 0:2n 
            distgivenk = Binomial(R, (1-ϵ)*k/(2n)+ϵ*(2n-k)/(2n))
            logprobs[k+1] = logpdf(distgivenk, d) + logpdf(distk, k)
        end 
        logpdftheoretical += counts[index] * logsumexp(logprobs)
    end 
    
    return logpdffunction, logpdftheoretical
end 

function testexpspace(n, ϵ, R, expsize)
    Random.seed!(1987)
    y = rand(0.1:0.1:0.9, expsize)
    q = Matrix{Float64}(undef, 2n+1, R+1)
    Random.seed!(1987)
    reads = rand(0:R, expsize)
    cm = countmap(zip(y, reads))
    y = [key[1] for key in keys(cm)]
    reads = [key[2] for key in keys(cm)]
    counts = collect(values(cm))
    Q, Qꜜ, hvals, binomialcoefficients = preparedata(n, y)
    
    logpdffunction = logprobderivedreadstestexpspace!(q, reads, R, R, counts, 0.0, 0.0, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
    
    logpdftheoretical = 0.0
    for index in eachindex(reads, counts, y)
        d = reads[index]
        distk = Binomial(2n, y[index])
        prob = 0.0 
        for k in 0:2n 
            distgivenk = Binomial(R, (1-ϵ)*k/(2n)+ϵ*(2n-k)/(2n))
            prob += pdf(distgivenk, d)*pdf(distk, k)
        end 
        logpdftheoretical += counts[index] * log(prob)
    end 
    
    return logpdffunction, logpdftheoretical
end 

function testcomparelikexpspacewiththeoritcallogspace(n, ϵ, R, expsize)
    Random.seed!(1987)
    y = rand(0.1:0.1:0.9, expsize)
    q = Matrix{Float64}(undef, 2n+1, R+1)
    Random.seed!(1987)
    reads = rand(0:R, expsize)
    cm = countmap(zip(y, reads))
    y = [key[1] for key in keys(cm)]
    reads = [key[2] for key in keys(cm)]
    counts = collect(values(cm))
    Q, Qꜜ, hvals, binomialcoefficients = preparedata(n, y)
    
    logpdffunction = logprobderivedreadstestexpspace!(q, reads, R, R, counts, 0.0, 0.0, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
    
    logpdftheoretical = 0.0
    logprobs = Vector{Float64}(undef, 2n+1)
    for index in eachindex(reads, counts, y)
        d = reads[index]
        distk = Binomial(2n, y[index])
        for k in 0:2n 
            distgivenk = Binomial(R, (1-ϵ)*k/(2n)+ϵ*(2n-k)/(2n))
            logprobs[k+1] = logpdf(distgivenk, d) + logpdf(distk, k)
        end 
        logpdftheoretical += counts[index] * logsumexp(logprobs)
    end 

    return logpdffunction, logpdftheoretical
end 


themax = 0.0
expsize = 100_000_000
for n in 1:10
    for R in 1:100
        for ϵ in 0.0:0.01:0.49
            ans1, ans2 = testlogspace(n, ϵ, R, expsize)
            m = abs(ans1-ans2)
            if m > themax 
                themax = m 
                #@show themax 
            end 
        end 
        #@show n R 
    end 
end # max voor werken in logspace: 6.571212907147128e-9 
themaxlogspace = themax  

themax = 0.0
for n in 1:10
    for R in 1:100
        for ϵ in 0.0:0.01:0.49
            ans1, ans2 = testexpspace(n, ϵ, R, expsize)
            m = abs(ans1-ans2)
            if m > themax 
                themax = m 
                #@show themax 
            end 
        end 
        #@show n R 
    end 
end # max voor werken in expspace: 5.664105628966354e-9. 
themaxexpspace = themax 

themax = 0.0
for n in 1:10
    for R in 1:100
        for ϵ in 0.0:0.01:0.49
            ans1, ans2 = testcomparelikexpspacewiththeoritcallogspace(n, ϵ, R, expsize)
            m = abs(ans1-ans2)
            if m > themax 
                themax = m 
                #@show themax 
            end 
        end 
        #@show n R 
    end 
end # max voor werken in expspace: 5.663878255290911e-9. 
themaxcompare = themax 

@show themaxlogspace themaxexpspace themaxcompare;

#=
Results with 1000_000
themaxlogspace = 6.782473064959049e-6
themaxexpspace = 5.754642188549042e-6
themaxcompare = 5.754642188549042e-6
=# 

##
## Test speed difference. 
##

for R in 1:100
    timedurationscalar = Millisecond(0.0)
    timedurationvector = Millisecond(0.0)
    for n in 1:10
        for ϵ in 0.1:0.1:0.4
            expsize = 1000_000 
            Random.seed!(1987)
            y = rand(0.0001:0.0001:0.9999, expsize)
            q = Matrix{Float64}(undef, 2n+1, R+1)
            Random.seed!(1987)
            reads = rand(0:R, expsize)
            cm = countmap(zip(y, reads))
            y = [key[1] for key in keys(cm)]
            reads = [key[2] for key in keys(cm)]
            counts = collect(values(cm))
            Rvec = fill(R, length(reads)) 
            Q, Qꜜ, hvals, binomialcoefficients = preparedata(n, y)
            before = now()
            for τC in 0.0:0.01:1.0 
                for τA in 0.0:0.01:1.0
                    logprobderivedreads!(q, reads, R, R, counts, τC, τA, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
                end 
            end 
            after = now()
            timedurationscalar += after - before

            uniquereads = unique(reads)
            q = Dict((d, R) => Vector{Float64}(undef, 2n+1) for d in uniquereads)
            before = now()
            for τC in 0.0:0.01:1.0 
                for τA in 0.0:0.01:1.0
                    logprobderivedreads!(q, reads, Rvec, [R], counts, τC, τA, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
                end 
            end 
            after = now()
            timedurationvector += after - before
        end # ϵ  
    end # n
    println(R, " takes ", timedurationscalar, " with scalar version, and ")
    println(R, " takes ", timedurationvector, " with vector version.")
    println()
    println()
end # R
