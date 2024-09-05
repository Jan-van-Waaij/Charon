using Primes 


let 
    prime = 2 
    while prime ≤ 10^9
        possibleprimes = [x^2-x+prime for x in 1:prime-1]
        if all(isprime.(possibleprimes))
            println(prime)
        end 
        prime = nextprime(prime+1)
    end 
end 

#=
N = 10000
primessmallerthanN = primes(N)

primesforwhichitworks = Int[] 
for prime in primessmallerthanN
    possibleprimes = [x^2-x+prime for x in 1:prime-1]
    if all(isprime.(possibleprimes))
        push!(primesforwhichitworks, prime)
    end 
end 

for prime in primesforwhichitworks
    println(prime)
end 
=#
