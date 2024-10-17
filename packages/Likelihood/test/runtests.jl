using Test, Likelihood, Distributions, GZip, LinearAlgebra, StatsBase, Random

using Likelihood: EigenExpansion, makeQ, makeQꜜ, calcmatrix

@testset "Likelihood" begin
    @testset "EigenExpansion and exp" begin
        x = rand(2)
        y = rand(3)
        A = Tridiagonal(x, y, x)
        F = EigenExpansion(A)
        @test A ≈ F.P*F.D*F.Pinv
        @test exp(Matrix(A)) ≈ exp(F,1.0)
        @test typeof(exp(F, 1.0)) == typeof(F.P)

        A = rand(2,2)
        B = A*A'
        F = EigenExpansion(B)
        @test typeof(exp(F,2.0)) == typeof(F.P)
        @test B ≈ F.P*F.D*F.Pinv
        @test exp(B) ≈ exp(F, 1.0)
        @test exp(2.0*B) ≈ exp(F, 2.0)
    end

    @testset "==" begin
        A = rand(2,2)
        B = A*A'
        C = copy(B)
        FB = EigenExpansion(B)
        FC = EigenExpansion(C)
        @test FB == FC 
        x = rand(3)
        y = rand(4)
        T = Tridiagonal(x, y, x)
        U = deepcopy(T)
        MT = Matrix(T)
        MU = Matrix(U)
        FT = EigenExpansion(T)
        FU = EigenExpansion(U)
        FMT = EigenExpansion(MT)
        FMU = EigenExpansion(MU)
        @test FT == FU
        @test FMT == FMU
    end

    @testset "makeQ" begin
        @test makeQ(1) == [0 1 0
                           0 -1 0
                           0 1 0]
    end
    
    @testset "makeQꜜ" begin
        @test makeQꜜ(1) == [0 3 0
                            0 -2 1
                            0 1 -2]
    end
    
    @testset "calcmatrix" begin
        for n in 1:10
            Q = makeQ(n)
            Qꜜ = makeQꜜ(n)
            FQ = EigenExpansion(Q)
            FQꜜ = EigenExpansion(Qꜜ)
            @test calcmatrix(FQ, FQꜜ, 0.0, 0.0) ≈ I
            τC = rand()
            τA = rand()
            M = calcmatrix(FQ, FQꜜ, τC, τA)
            @test typeof(FQ.P) == typeof(M)
            @test M == exp(FQ, τA)*exp(FQꜜ, τC)
            M = calcmatrix(FQ, FQꜜ, 0.0, τA)
            @test M ≈ exp(FQ, τA)
            M = calcmatrix(FQ, FQꜜ, τC, 0.0)
            @test M ≈ exp(FQꜜ, τC)
        end 
    end

    @testset "preparedata" begin
        y = 0.1:0.1:0.9
        for n in 1:10
            t = preparedata(n, y)
            @test length(t) == 4 
            @test t isa NamedTuple 
            Q, Qꜜ, hvals, binomialcoefficients = t
            @test Q == EigenExpansion(makeQ(n))
            @test Qꜜ == EigenExpansion(makeQꜜ(n))
            @test size(hvals) == (2n+1, length(y))
            for k in 0:2n 
                for yindex in eachindex(y) 
                    @test hvals[k+1, yindex] == y[yindex]^k*(1-y[yindex])^(2n-k)
                end 
            end 
            @test binomialcoefficients == binomial.(2n, 0:2n)
        end 
    end

    @testset "logprobderivedreads!" begin
        # when τC=τA=0, then logprobderivedreads! calculates the log probability 
        # of the binomial distribution with parameters R and (1-ϵ)*k/2n+ϵ*k/(2n), where k itself is binomially distributed with parameters 2n and y. 
        for n in 1:10
            for R in 1:100
                for ϵ in 0.1:0.1:0.4
                    expsize = 1000 
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
                    τC = 0.3 
                    τA = 0.2 
                    @test logprobderivedreads!(q, reads, R, R, counts, -0.1, τA, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals) == -Inf 
                    @test logprobderivedreads!(q, reads, R, R, counts, τC, -0.1, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals) == -Inf 
                    @test logprobderivedreads!(q, reads, R, R, counts, τC, τA, -0.1, n, Q, Qꜜ, binomialcoefficients, hvals) == -Inf 
                    @test logprobderivedreads!(q, reads, R, R, counts, τC, τA, 1.1, n, Q, Qꜜ, binomialcoefficients, hvals) == -Inf 
                    logpdffunction = logprobderivedreads!(q, reads, R, R, counts, 0.0, 0.0, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
                    #@assert logpdffunction > -Inf 
                    logpdftheoretical = 0.0
                    for index in eachindex(reads, counts, y)
                        d = reads[index]
                        distk = Binomial(2n, y[index])
                        prob = 0.0 
                        for k in 0:2n 
                            distgivenk = Binomial(R, (1-ϵ)*k/(2n)+ϵ*(2n-k)/(2n))
                            prob += pdf(distgivenk, d)*pdf(distk, k)
                        end 
                        @assert prob > 0.0
                        @assert prob ≤ 1.0 
                        logpdftheoretical += counts[index] * log(prob)
                    end 
                    @test logpdffunction ≈ logpdftheoretical 
                    println("Difference: ", abs(logpdffunction - logpdftheoretical))
                end # ϵ  
            end # R
        end # n 

        # logprobderivedreads! with R or fill(R, length(reads)) should be the same. 
        for n in 1:10
            for R in 1:100 
                ϵ = rand(Uniform(0.0, 0.5))
                y = rand(1000)
                Rvec = fill(R, 1000)
                uniqueRvec = [R]
                reads = rand(0:R, 1000)
                cm = countmap(zip(reads, y))
                reads = [key[1] for key in keys(cm)]
                y = [key[2] for key in keys(cm)]
                counts = collect(values(cm))
                Rvec = fill(R, length(counts))
                Q, Qꜜ, hvals, binomialcoefficients = preparedata(n, y)
                qmatrix = Matrix{Float64}(undef, 2n+1, R+1)
                qdic = Dict((d, R) => Vector{Float64}(undef, 2n+1) for d in unique(reads))
                for τC in 0.0:0.1:1.0
                    for τA in 0.0:0.1:1.0
                        for ϵ in 0.1:0.1:0.4
                            resultintversion = logprobderivedreads!(qmatrix, reads, R, R, counts, τC, τA, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
                            resultvecversion = logprobderivedreads!(qdic, reads, Rvec, uniqueRvec, counts, τC, τA, ϵ, n, Q, Qꜜ, binomialcoefficients, hvals)
                            @test resultintversion == resultvecversion
                        end # ϵ
                    end #τA
                end #τC
            end # R 
        end # n 
    end # testset logprobderivedreads!
end # testset Likelihood