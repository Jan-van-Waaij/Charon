using Likelihood: calcmatrix, makeQ, makeQꜜ, exp, EigenExpansion

n = 10

Q = EigenExpansion(makeQ(n))
Qꜜ = EigenExpansion(makeQꜜ(n))

τrange = 0.0001:0.0001:1.0

@time for τA in τrange
    expτAQ = exp(Q, τA)
    for τC in τrange
        expτCQꜜ = exp(Qꜜ, τC)
        M = expτAQ * expτCQꜜ
    end
end