
println(length(ARGS)) 
println(typeof(ARGS[1]))

n = parse()




#= 
function Q(n)
    lowerdiagonal = [0.5*i*(i-1) for i in 1:n]
    diagonal = Float64[-i*(n-i) for i in 0:n]
    upperdiagonal = [0.5*(n-i)*(n-i-1) for i in 0:n-1]
    return diagm(-1 => lowerdiagonal, 0 => diagonal, 1=> upperdiagonal)
end

function Qꜜ(n)
    lowerdiagonal = [0.5*i*(i-1) for i in 1:n]
    diagonal = Float64[-i*(n-i+1) for i in 0:n]
    upperdiagonal = [0.5*(n-i+1)*(n-i) for i in 0:n-1]
    return diagm(-1 => lowerdiagonal, 0 => diagonal, 1=> upperdiagonal)
end 

generatehartikel(y, n) = [y^k*(1-y)^(n-k) for k in 0:n]

generatehvergissing(y, n) = [binomial(n, k)*y^k*(1-y)^(n-k) for k in 0:n]

partikel(t1, t2, n, y) = exp(Q(n)*t2)*exp(Qꜜ(n)*t1)*generatehartikel(y, n)

pvergissing(t1, t2, n, y) = exp(Q(n)*t2)*exp(Qꜜ(n)*t1)*generatehvergissing(y, n)

for n in 1:10
    for y in 0.1:0.1:1.0
        for t1 in 0.1:0.1:1.0
            for t2 in 0.1:0.1:1.0
                @test all(partikel(t1, t2, n, y) .  )
            end 
        end 
    end 
end 

=#
