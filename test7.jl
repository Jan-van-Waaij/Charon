using Likelihood, Distributions, LinearAlgebra

function probderivedalleles!(p::AbstractVector{Float64}, h::AbstractVector{Float64}, M::AbstractMatrix{Float64}, binomialcoefficients::AbstractVector{Int})
    # mul! is 60% faster. 
    mul!(p, M, h) # equivalent to p = M*h, but in place, without allocating new vectors.
    p .*= binomialcoefficients # in place computation of p = p*binomialcoefficients
    p .= max.(0.0, p)
    return p  
end

for n in 11:30
    println("starting with ", n)
    Q = EigenExpansion(Matrix(Likelihood.makeQ(n)))
    Qꜜ = EigenExpansion(Matrix(Likelihood.makeQꜜ(n)))
    p = Vector{Float64}(undef, 2n+1)
    binomialcoefficients = binomial.(2n, 0:2n)
    sumdiffmax = 0.0
    pmin = 0.0 
    for τA in 0.01:0.01:1.0
        for τC in 0.01:0.01:1.0
            M = Likelihood.calcmatrix(Q, Qꜜ, τC, τA)
            for y in 0.01:0.01:1.0-0.01
                h = [y^k*(1-y)^(2n-k) for k in 0:2n]
                probderivedalleles!(p, h, M, binomialcoefficients)
                diff = abs(sum(p) - 1.0)
                mi = minimum(p)
                if diff > sumdiffmax 
                    sumdiffmax = diff 
                end 
                if mi < pmin 
                    pmin = mi 
                end 
            end 
        end 
    end 
    @show sumdiffmax pmin 
    println("Finished with ", n)
end 



