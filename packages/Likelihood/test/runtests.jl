using Test, Likelihood, Distributions, GZip, LinearAlgebra

using Likelihood: EigenExpansion, makeQ, makeQꜜ, calcmatrix, probderivedalleles!

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
        M = Matrix(T)
        FT = EigenExpansion(T)
        FU = EigenExpansion(U)
        FM = EigenExpansion(M)
        @test FT ≈ FM # SMatrix can make use of different methods
        # for eigen decomposition. 
        @test FT == FU 
    end
    
    @testset "≈" begin
        A = rand(2,2)
        B = A*A'
        C = copy(B)
        C[1,1] += eps(C[1,1])
        @test B != C
        FB = EigenExpansion(B)
        FC = EigenExpansion(C)  
        @test FB ≈ FC 
        x = rand(3)
        y = rand(4)
        z = copy(y)
        z[1] += eps(z[1])
        T = Tridiagonal(x, y, x)
        U = Tridiagonal(x, z, x)
        FT = EigenExpansion(T)
        FU = EigenExpansion(U)
        @test FT ≈ FU
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
        Q = makeQ(2)
        Qꜜ = makeQꜜ(2)
        FQ = EigenExpansion(Q)
        FQꜜ = EigenExpansion(Qꜜ)
        τC = rand()
        τA = rand()
        M = calcmatrix(FQ, FQꜜ, τC, τA)
        @test M == exp(FQ, τA)*exp(FQꜜ, τC)
        @test calcmatrix(FQ, FQꜜ, 0.0, 0.0) ≈ I
        @test typeof(FQ.P) == typeof(M)
    end
    @testset "probderivedalleles!" begin
        for n in 1:10
            @show n
            p = Vector{Float64}(undef, 2n+1)
            for y in 0.1:0.1:0.9
                h = [y^k*(1.0-y)^(2n-k) for k in 0:2n]
                for τC in 0.1:0.1:10.0
                    for τA in 0.1:0.1:10.0
                        Q = makeQ(n)
                        Qꜜ = makeQꜜ(n)
                        FQ = EigenExpansion(Q)
                        FQꜜ = EigenExpansion(Qꜜ)
                        M = calcmatrix(FQ, FQꜜ, τC, τA)
                        binomialcoefficients = binomial.(2n, 0:2n)
                        q = probderivedalleles!(p, h, M, binomialcoefficients)
                        @test p == q
                        @test isprobvec(p)
                        @test p == (M*h).*binomialcoefficients 
                    end 
                end
            end 
        end
    end
    
    @testset "preparedatafixedn" begin
        for n in 1:10
            R = 20
            y = rand(50)
            Q, Qꜜ, hvals, binomialcoefficients, p, binomialprobs = preparedatafixedn(n, R, y)
            @test Q == EigenExpansion(makeQ(n))
            @test Qꜜ == EigenExpansion(makeQꜜ(n))
            for i in eachindex(y)
                for k in 0:2n 
                    @test hvals[k+1, i] == y[i]^k*(1.0-y[i])^(2n-k)
                end 
            end 
            @test binomialcoefficients == binomial.(2n, 0:2n)
            @test length(p) == 2n+1 
        end
        for n in 1:10
            for R in [1, 2, 3, 10, 30, 100]
                for d in 0:R
                    y = rand(50)
                    Q, Qꜜ, hvals, binomialcoefficients, p, binomialprobs = preparedatafixedn(n, R, y)
                    for k in 0:2n 
                        @test binomialprobs[k+1, d+1] == pdf(Binomial(R, k/2n), d)
                    end
                end  
            end
        end
    end

    @testset "logprobderivedreads!" begin
        n = 4
        R = 20
        y = rand(50)
        Q, Qꜜ, hvals, binomialcoefficients, p, binomialprobs = preparedatafixedn(n, R, y)
        d = rand(0:R, length(y))
        τC = 0.0
        τA = 0.1 
        @test logprobderivedreads!(d, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals) == -Inf
        τC = 0.1
        τA = 0.0 
        @test logprobderivedreads!(d, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals) == -Inf
        τC = 0.0
        τA = 0.0 
        @test logprobderivedreads!(d, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals) == -Inf
        τC = -0.1
        τA = 0.1 
        @test logprobderivedreads!(d, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals) == -Inf
        τC = 0.1
        τA = -0.1 
        @test logprobderivedreads!(d, τC, τA, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals) == -Inf
        n = 4
        R = 2
        y = rand(3)
        ds =  [[a, b, c] for a in 0:2 for b in 0:2 for c in 0:2]
        logprobs = similar(ds, Float64)
        Q, Qꜜ, hvals, binomialcoefficients, p, binomialprobs = preparedatafixedn(n, R, y)
        for i in eachindex(ds, logprobs)
            logprobs[i] = logprobderivedreads!(ds[i], 0.2, 0.1, Q, Qꜜ, p, binomialcoefficients, binomialprobs, hvals)
        end 
        probs = exp.(logprobs)
        @test isprobvec(probs) 
    end
end