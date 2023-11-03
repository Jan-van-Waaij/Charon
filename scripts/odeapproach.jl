using Likelihood: makeQ, makeQꜜ
using LinearAlgebra, Distributions, SpecialFunctions

#=
Here we take advantage of the fact that x(t)=e^{At}y is a solution the the ordinary stochastic differential equation
dx(t) = At
x(0) = y. 

We want to test the speed of this approach, how close it is to the usual approach, it's speed, whether it results in a probability vector and if the values are positive.  
=# 

function probderivedreadstandard(n::Int, τC::Float64, τA::Float64, y::Float64)
    Q = Matrix(makeQ(n))
    Qꜜ = Matrix(makeQꜜ(n))
    p = exp(Q*τA)*exp(Qꜜ*τC)*[y^k*(1-y)^(2n-k) for k in 0:2n]
    for k in 0:2n
        p[k+1] = binomial(2n,k)*p[k+1]
    end 
    return p 
end 

function probderivedreadsodeapproach(n::Int, τC::Float64, τA::Float64, y::Float64, nsteps::Int=1000_000)
    Qꜜ = makeQꜜ(n)
    Q = makeQ(n)
    p = [y^k*(1-y)^(2n-k) for k in 0:2n]
    h = τC/nsteps
    k1 = similar(p)
    k2 = similar(p)
    k3 = similar(p)
    k4 = similar(p)
    for k in Base.OneTo(nsteps)
        mul!(k1, Qꜜ, p)
        mul!(k2, Qꜜ, p + (h/2).*k1)
        mul!(k3, Qꜜ, p + (h/2).*k2)
        mul!(k4, Qꜜ, p + h.*k3)
        p .= p .+ (h/6).*(k1.+k2.+2k3.+k4)
    end 
    h = τA/nsteps
    for k in Base.OneTo(nsteps)
        mul!(k1, Q, p)
        mul!(k2, Q, p .+ (h/2).*k1)
        mul!(k3, Q, p .+ (h/2).*k2)
        mul!(k4, Q, p .+ h.*k3)
        p .= p .+ (h/6).*(k1.+k2.+2k3.+k4)
    end
    for k in eachindex(p)
        p[k] = exp(logabsbinomial(2n,k-1)[1])*p[k]
    end 
    return p 
end 

τrange = 0.001:0.001:2.0
yrange = 0.001:0.001:0.999
for n in 1:200
    for _ in Base.OneTo(100)
        τC = rand(τrange)
        for _ in Base.OneTo(100)
            τA = rand(τrange)
            for _ in Base.OneTo(100)
                y = rand(yrange)
                p = probderivedreadsodeapproach(n, τC, τA, y, 100_000)
                if !isprobvec(p)
                    @show n τC τA y 
                    error()
                end 
            end
        end
    end 
    println("finished with n = ", n)
end 