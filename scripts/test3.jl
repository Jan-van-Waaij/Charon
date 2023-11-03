using Likelihood, LinearAlgebra, Distributions
using Likelihood: calcmatrix, makeQ, makeQꜜ

function test(p, M, y)   
    n = (size(M,1)-1)÷2     
    h = [y^k*(1-y)^(2n-k) for k in 0:2n]
    binomialcoefficients = binomial.(2n, 0:2n)
    mul!(p, M, h)
    p .*= binomialcoefficients
    p .= max.(0.0, p)
    return p 
end

function overtest(Q, Qꜜ, p, τCrange, τArange, yrange)
    n = (size(Q.P,1)-1)÷2
    for τC in τCrange
        for τA in τArange
            M = calcmatrix(Q, Qꜜ, τC, τA)
            for y in yrange 
                test(p, M, y)
                @assert all(p .≥ 0.0)
                if !isprobvec(p)
                    @show τC, τA, y, sum(p)
                    #error("oeps")
                end
            end 
        end 
    end 
end 

function overovertest(n, range)
    println("testing ", n)
    Q = EigenExpansion(makeQ(n))
    Qꜜ = EigenExpansion(makeQꜜ(n))
    p = Vector{Float64}(undef,2n+1)
    overtest(Q, Qꜜ, p, range, range, range)
end 

nrange = 14:28
range = 0.1:0.1:0.9

for n in nrange
    overovertest(n, range)
end 

println("finished")