using Likelihood, LinearAlgebra, Distributions, StaticArrays
using Likelihood: calcmatrix, makeQ, makeQꜜ
function sloppytesteigenexpansion()
    τrange = 1e-6:1e-6:2.0
    yrange = 1e-6:1e-6:1.0-1e-6
    for n in 1:100
        Q = EigenExpansion(makeQ(n))
        Qꜜ = EigenExpansion(makeQꜜ(n))
        binomialcoefficients = SVector{2n+1,Int}(binomial.(2n,0:2n))
        p = MVector{2n+1,Float64}(undef) 
        for _ in Base.OneTo(100)
            τC = rand(τrange)
            τA = rand(τrange)
            M = calcmatrix(Q, Qꜜ, τC, τA)
            for _ in Base.OneTo(100)
                y = rand(yrange)
                h = SVector{2n+1,Float64}([y^k*(1-y)^(2n-k) for k in 0:2n])
                mul!(p, M, h)
                p .*= binomialcoefficients
                p .= max.(p,0.0)
                if !(1-1e-6 < sum(p) < 1+1e-6)
                    @show n τC τA y 
                    throw(AssertionError("p is not a probability vector"))
                end 
            end 
        end 
        println("EigenExpansion test, finished with ", n)
    end # gaat goed tot en met n ≤ 13
end 

sloppytesteigenexpansion()

#=
Fails with 
n = 16
τC = 1.923546
τA = 0.620386
y = 0.931598
=#

