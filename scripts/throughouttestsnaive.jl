using Likelihood, LinearAlgebra, Distributions, StaticArrays, SpecialFunctions
using Likelihood: makeQ, makeQꜜ

function throughouttestnaive()
    for n in 1:70
        Q = SizedMatrix{2n+1,2n+1,Float64}(Matrix(makeQ(n)))
        Qꜜ = SizedMatrix{2n+1,2n+1,Float64}(Matrix(makeQꜜ(n)))
        binomialcoefficients = SVector{2n+1,Float64}(exp.(logabsbinomial.(2n, k)[1] for k in 0:2n))
        τrange = 0.0:0.001:2.0
        yrange = 0.0:0.001:1.0
        for y in yrange
            h = SVector{2n+1,Float64}([y^k*(1-y)^(2n-k) for k in 0:2n])
            for τC in τrange
                for τA in τrange
                    p = exp(τC*Qꜜ)*h
                    p = exp(τA*Q)*p 
                    p = binomialcoefficients .* p
                    if !isprobvec(p)
                        @show n y τC τA 
                        throw(AssertionError("p is not a probability vector"))
                    end
                end 
            end 
        end 
        println("finished with ", n)
    end 
end 

throughouttestnaive()

#=
Fails with 
n = 7
y = 1.0
τC = 0.0
τA = 0.027
=#