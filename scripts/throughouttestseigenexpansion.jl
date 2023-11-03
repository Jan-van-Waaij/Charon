using Likelihood, LinearAlgebra, Distributions, StaticArrays
using Likelihood: calcmatrix, makeQ, makeQꜜ
function throughouttesteigenexpansion()
    τrange = 0.0:0.001:2.0
    yrange = 0.0:0.001:1.0
    for n in 1:13
        Q = EigenExpansion(makeQ(n))
        Qꜜ = EigenExpansion(makeQꜜ(n))
        binomialcoefficients = SVector{2n+1,Int}(binomial.(2n,0:2n))
        for τC in τrange
            for τA in τrange
                M = calcmatrix(Q, Qꜜ, τC, τA)
                p = MVector{2n+1,Float64}(undef) 
                for y in yrange
                    h = SVector{2n+1,Float64}([y^k*(1-y)^(2n-k) for k in 0:2n])
                    mul!(p, M, h)
                    p .*= binomialcoefficients
                    p .= max.(p,0.0)
                    if !isprobvec(p)
                        @show n τC τA y 
                        throw(AssertionError("p is not a probability vector"))
                    end 
                end 
            end 
        end 
        println("EigenExpansion test, finished with ", n)
    end # gaat goed tot en met n ≤ 13
end 

throughouttesteigenexpansion()
# Fails for n = 13, τC = τA = 0.0 and y=0.984. 

# works for n ≤ 13 with y 1.0e-6:1.0e-6:(1-1.0e-6) 
# en τ 0.01:0.01:2.0.
