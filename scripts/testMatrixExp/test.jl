using LinearAlgebra
using Expokit 
using Dates 
using Likelihood: makeQ, makeQꜜ 

nrange = 1:10
numexp = 100_000 
svec = rand(numexp)


for n in nrange
    @show n  
    starttime = now()
    M = makeQ(n)
    Mꜜ = makeQꜜ(n)
    F = EigenExpansion(M)
    Fꜜ = EigenExpansion(Mꜜ)
    for M in [F, Fꜜ]
        for index in Base.OneTo(numexp)
            v = rand(2n+1)
            exp(M, svec[index])*v 
        end 
    end 
    finishtime = now()
    println("Time my method: ", finishtime-starttime)
    
    starttime = now()
    for MT in [M, Mꜜ]
        for index in Base.OneTo(numexp)
            v = rand(2n+1)
            expmv!(svec[index], MT, v)
        end 
    end 
    finishtime = now()
    println("Time my method: ", finishtime-starttime)
end 

for n in nrange
    @show n  
    M = makeQ(n)
    Mꜜ = makeQꜜ(n)
    F = EigenExpansion(M)
    Fꜜ = EigenExpansion(Mꜜ)
    aboutthesame = BitVector(undef, numexp)
    for index in eachindex(Base.OneTo(numexp), aboutthesame)
        v = rand(2n+1)
        ansmymethod = exp(F, svec[index])*v 
        ansexpokitmethod = expmv(svec[index], M, v)
        aboutthesame[index] = ansmymethod ≈ ansexpokitmethod
    end 
    println(length(aboutthesame)-sum(aboutthesame))
    println(mean(aboutthesame))
    for index in eachindex(Base.OneTo(numexp), aboutthesame)
        v = rand(2n+1)
        ansmymethod = exp(Fꜜ, svec[index])*v 
        ansexpokitmethod = expmv(svec[index], Mꜜ, v)
        aboutthesame[index] = ansmymethod ≈ ansexpokitmethod
    end 
    println(length(aboutthesame)-sum(aboutthesame))
    println(mean(aboutthesame))
end 

for n in nrange
    @show n  
    M = makeQ(n)
    M_Matrixtype = Matrix(M)
    Mꜜ = makeQꜜ(n)
    Mꜜ_Matrixtype = Matrix(Mꜜ)
    F = EigenExpansion(M)
    Fꜜ = EigenExpansion(Mꜜ)
    aboutthesame = BitVector(undef, numexp)
    for index in eachindex(Base.OneTo(numexp), aboutthesame)
        v = rand(2n+1)
        #ansmymethod = exp(F, svec[index])*v
        ansstandardmethod = exp(svec[index]*M_Matrixtype)*v  
        ansexpokitmethod = expmv(svec[index], M, v)
        aboutthesame[index] = ansexpokitmethod ≈ ansstandardmethod
    end 
    println(length(aboutthesame)-sum(aboutthesame))
    println(mean(aboutthesame))
    for index in eachindex(Base.OneTo(numexp), aboutthesame)
        v = rand(2n+1)
        #ansmymethod = exp(Fꜜ, svec[index])*v
        ansstandardmethod = exp(svec[index]*Mꜜ_Matrixtype)*v  
        ansexpokitmethod = expmv(svec[index], M, v)
        aboutthesame[index] = ansexpokitmethod ≈ ansstandardmethod
    end 
    println(length(aboutthesame)-sum(aboutthesame))
    println(mean(aboutthesame))
end 

istriuvec = BitVector(undef, numexp)
for n in nrange 
    for index in eachindex(isdiagvec)
        M = diagm(0 => rand(2n+1), 1 => rand(2n))
        E = exp(M)
        istriuvec[index] = istriu(E)
    end
    println(length(istriuvec) - sum(istriuvec))
    println(mean(istriuvec))
end  

function makeorthogonalmatrix(n::Int, M::Matrix{Float64})
    n ≥ 1 || error("n should be one or larger.")
    size(M) == (2,2) || error("M should be 2x2 matrix")
    BM = zeros(2n,2n)
    for i in 0:n-1 
        BM[2i+1:2i+2, 2i+1:2i+2] .= M  
    end 
    return BM 
end 

M = [0.6   0.8
0.8  -0.6]

n = 10

BM = makeorthogonalmatrix(n, M)
X = BM' * diagm(1:2n) * BM 

TRUEEXPX = BM' * diagm(exp.(1:2n)) * BM 

v = rand(2n)
exp(X)*v ≈ TRUEEXPX*v

F = EigenExpansion(X)
myexpsol = exp(F, 1.0)

myexpsol*v ≈ TRUEEXPX*v

expmv(1.0, X, v) ≈ TRUEEXPX*v 