# runmcmconallfiles.jl
# Author: Jan van Waaij
# Created: 4 May 2023
# Last update: 4 May 2023
#
# Purpose:
# Run MCMC.
#
# Installation:
# Please read installinstructions.txt.
# 
# Usage:
# open file as julia path/to/runmcmconallfiles.jl
#
# Requirement:
# Frequency files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/frequencies/freqgenotypes
# Read files in /home/projects/SedaDNA_Demo/simulationsMSPRIME/data/reads
# 
# Output:
# ????
# The CSV file has 10 columns: NeAnchorPop, NeAncientPop, jointtimeyears, trueτC, trueτA, truenumind, τC, τA, numind, loglik.

#using Pkg
#Pkg.update()

using Data, # For parsing frequency and read files. 
      Distributions, # To specify the prior.
      GZip, # for opening zip files.
      Likelihood # the module we are testing

vcffilespath = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "vcffiles")
allfiles = readdir(vcffilespath)[2:end] #["two_demes_3000_6000_30000_0.vcf.gz"]



fileswithoutextensionvec = [file[begin:end-7] for file in allfiles]

numberofindividualsinthesedimentvec =  [1, 2, 4] 
nmax = 4

firsttimeold = true
firsttimenew = true
firsttimeneweigenexpansion = true 
firsttimesupernew = true 

Qs = [EigenExpansion(makeQ(k)) for k in 1:nmax]
# A vector of Qꜜ matrices. 
Qꜜs = [EigenExpansion(makeQꜜ(k)) for k in 1:nmax]
# A vector with quantities y[j]^d*(1.0-y[j])^(2n-d).
# Binomial coefficients 2n choose k, k=0,...,2n, and then concatenated for all n=1,...,nmax.
binomialcoefficientsunknownn = vcat((binomial.(2n, 0:2n) for n in 1:nmax)...)
# The auxiliary vector for the function. 
punknownn = Vector{Float64}(undef, 2*nmax+1)
# binomialprobs is a nmax*(nmax+2)x(R+1)-matrix, where binomialprobs[i,j] contains the 
# probability of d successes out of R experiments, with probability of success equal to 
# k/2n, where d=0,...,R, k=0,...,2n and i = n^2 + k and j = d+1.  
binomialprobsunknownn = Matrix{Float64}(undef, nmax*(nmax+2), 100+1)
for n in 1:nmax
    binomialdistributions = Binomial.(100, range(0.0, 1.0, length = 2n+1))
    for d in 0:100    
        for k in 0:2n
            i = n^2+k 
            j = d+1 
            binomialprobsunknownn[i,j] = pdf(binomialdistributions[k+1], d)
        end 
    end 
end


for filewithoutextension in fileswithoutextensionvec
    # filewithoutextension is of the form two_demes_Ne1_Ne2_jointimeyears_0
    frequencyfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "frequencies", 
    "freqanchor", string(filewithoutextension, "_freqanchor.txt.gz"))
    gzopen(frequencyfile) do f 
        global freqs = parsefreqfile(f)
    end
    acceptedindices = 0.0 .< freqs .< 1.0
    global freqs = freqs[acceptedindices]
    hvalsunknownn = Matrix{Float64}(undef, nmax*(nmax+2), length(freqs))
    for j in eachindex(freqs)
        for n in 1:nmax 
            for d in 0:2n 
                i = (n-1)*(n+1) + d + 1
                hvalsunknownn[i, j] = freqs[j]^d*(1.0-freqs[j])^(2n-d)
            end 
        end 
    end 

    for truenumind in numberofindividualsinthesedimentvec
        readsfile = joinpath("/", "home", "projects", "SedaDNA_Demo", "simulationsMSPRIME", "data", "reads", 
            string(filewithoutextension, '_', truenumind, "_100_reads.txt.gz") )
        if !isfile(readsfile)
            println(readsfile, " does not exist.")
            continue
        end 
        gzopen(readsfile) do f 
            global derivedreadsvec, coveragevec = parsereadfile(f)
        end
        global derivedreadsvec = derivedreadsvec[acceptedindices]
        global coveragevec = coveragevec[acceptedindices]
        @assert length(freqs) == length(derivedreadsvec) == length(coveragevec)
        @assert length(unique(coveragevec)) == 1
        coverage = coveragevec[1]

        println("==============================")
        println("new set of tests")
        @show filewithoutextension truenumind 
        #=
        println("old version")
        @time for τC in 0.1:0.001:1.0
            for τA in 0.1:0.01:1.0
                logprobderivedreads(coverage, derivedreadsvec, truenumind, τC, τA, freqs)
            end 
        end
        println("new version")
        =#
        n = truenumind
        R = coverage
        y = freqs
        Q = Matrix(makeQ(n))
        Qꜜ = Matrix(makeQꜜ(n))
        eigenQ = EigenExpansion(Q)
        eigenQꜜ = EigenExpansion(Qꜜ)
        p = Vector{Float64}(undef, 2n+1)
        h = Vector{Float64}(undef, 2n+1)
        binomialcoefficients = binomial.(2n, 0:2n)
        binomialdistributions = Binomial.(R, range(0.0, 1.0, length = 2n+1)) # with Distributions package. 
        binomialprobs = Matrix{Float64}(undef, 2n+1, R+1)
        for d in 0:R 
            for k in 1:2n+1
                binomialprobs[k, d+1] = pdf(binomialdistributions[k], d)
            end 
        end 
        hvals = Matrix{Float64}(undef, 2n+1, length(y))
        for j in eachindex(y)
            for k in 0:2n 
                hvals[k+1, j] = y[j]^k*(1.0-y[j])^(2n-k)
            end 
        end
        #=
        for τC in 0.1:0.01:1.0
            for τA in 0.1:0.01:1.0
                oldlik = logprobderivedreads(coverage, derivedreadsvec, truenumind, τC, τA, freqs)
                newlik = logprobderivedreads!(derivedreadsvec, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
                altlik = logprobderivedreads!(derivedreadsvec, τC, τA, Q, Qꜜ, p, h, binomialcoefficients, binomialprobs, freqs)
                @assert oldlik ≈ altlik
                @assert newlik ≈ altlik
            end 
        end
        @time for τC in 0.1:0.001:1.0
            for τA in 0.1:0.01:1.0
                logprobderivedreads!(derivedreadsvec, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
            end 
        end
        =#
        println("old version")
        @time for τC in 0.1:0.01:1.0
            for τA in 0.1:0.01:1.0
                logprobderivedreads(coverage, derivedreadsvec, truenumind, τC, τA, freqs)
            end 
        end
        if firsttimeold
            @code_warntype logprobderivedreads(coverage, derivedreadsvec, truenumind, 0.2, 0.3, freqs)
            global firsttimeold = false
        end 
        println("new version")
        @time for τC in 0.1:0.01:1.0
            for τA in 0.1:0.01:1.0
                logprobderivedreads!(derivedreadsvec, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
            end 
        end
        if firsttimenew
            @code_warntype logprobderivedreads!(derivedreadsvec, 0.2, 0.3, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
            global firsttimenew = false 
        end 
        println("new version with EigenExpansion")
        @time for τC in 0.1:0.01:1.0
            for τA in 0.1:0.01:1.0
                logprobderivedreads!(derivedreadsvec, τC, τA, eigenQ, eigenQꜜ, p, binomialcoefficients, binomialprobs, hvals)
            end 
        end
        if firsttimeneweigenexpansion
            @code_warntype logprobderivedreads!(derivedreadsvec, 0.2, 0.3, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
            global firsttimeneweigenexpansion = false 
        end 
        println("super new version")
        @time for τC in 0.1:0.01:1.0
            for τA in 0.1:0.01:1.0
                altlik = logprobderivedreads!(derivedreadsvec, n, τC, τA, Qs, Qꜜs, punknownn, binomialcoefficientsunknownn , binomialprobsunknownn, hvalsunknownn)
            end 
        end
        if firsttimesupernew
            @code_warntype logprobderivedreads!(derivedreadsvec, n, 0.2, 0.3, Qs, Qꜜs, punknownn, binomialcoefficientsunknownn , binomialprobsunknownn, hvalsunknownn)
            global firsttimesupernew = false 
        end 
        println("check correctness")
        @time for τC in 0.1:0.01:1.0
            for τA in 0.1:0.01:1.0
                oldlik = logprobderivedreads(coverage, derivedreadsvec, truenumind, τC, τA, freqs)
                newlik = logprobderivedreads!(derivedreadsvec, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
                supernewlik = logprobderivedreads!(derivedreadsvec, n, τC, τA, Qs, Qꜜs, punknownn, binomialcoefficientsunknownn , binomialprobsunknownn, hvalsunknownn)
                @assert oldlik ≈ newlik
                @assert oldlik ≈ supernewlik
                @assert newlik ≈ supernewlik
            end 
        end
    end # truenumind
end # for filewithoutextension