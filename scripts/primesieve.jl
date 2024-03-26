# On the cluster, Julia is installed at 
# /home/ctools/julia-1.6.7/bin/julia
# You can install Primes, Test and BenchmarkTools as follows:
#using Pkg
#Pkg.add(["Primes", "Test", "BenchmarkTools"])
using BenchmarkTools, Base.Threads, Primes, Test, BenchmarkTools, Base.Iterators
 # look at Iterators.partition
function primesieve(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    for i in 2:isqrt(n)
        if sieve[i] # i is prime 
            view(sieve, 2*i:i:length(sieve)) .= false # number 2*i, 3*i, ... are not prime. 
        end 
    end 
    return (1:n)[sieve]
end 

function primesievethreaded1(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    sqrtn = isqrt(n) 
    N = sqrtn ÷ (2nthreads()-2) 
    for i in 2:sqrtn
        if sieve[i] # i is prime
            pieces = partition(2*i:i:length(sieve), N)
            @threads for piece in collect(pieces) 
               view(sieve, piece) .= false # number 2*i, 3*i, ... are not prime.
            end  
        end 
    end 
    return (1:n)[sieve]
end 

function primesievethreaded2(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime.
    numthreads = nthreads() 
    sqrtn = isqrt(n) 
    for i in 2:sqrtn
        if sieve[i] # i is prime
            range = 2*i:i:length(sieve)
            N = length(range) ÷ (numthreads-1)
            pieces = partition(range, N)
            @threads for piece in collect(pieces) 
               view(sieve, piece) .= false # number 2*i, 3*i, ... are not prime.
            end  
        end 
    end 
    return (1:n)[sieve]
end 

function primesievethreaded2(n::Integer, k::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime.
    numthreads = nthreads() 
    sqrtn = isqrt(n) 
    for i in 2:sqrtn
        if sieve[i] # i is prime
            range = 2*i:i:length(sieve)
            N = length(range) ÷ k
            pieces = partition(range, N)
            @threads for piece in collect(pieces) 
               view(sieve, piece) .= false # number 2*i, 3*i, ... are not prime.
            end  
        end 
    end 
    return (1:n)[sieve]
end 

println("num threads ", nthreads())
trueprimes = primes(1000_000_000)
for m in 1:10*nthreads()
    for _ in 1:10
        println("k is ", m )
        @test primesievethreaded2(1000_000_000, m) == trueprimes
    #@btime primesievethreaded2(1000_000_000, m) seconds=600
    end 
end 

    
maxnum = 100_000
trueprimes = primes(maxnum)
@test primesieve(maxnum) == trueprimes
@test primesievethreaded1(maxnum) == trueprimes
@test primesievethreaded2(maxnum) == trueprimes

maxnum = 1000_000_000
@btime primesieve($maxnum) seconds=600;
@btime primesievethreaded1($maxnum) seconds=600;
@btime primesievethreaded2($maxnum) seconds=600;

# start julia with julia --threads=auto or julia --threads=14 # or another number. 
#=
function primesievemultithreaded(n::Integer)
    numthreads = nthreads() - 1 
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    for i in 2:isqrt(n)
        if sieve[i] # i is prime 
            # divide work.
            maxmultiple = n÷i 
            divide = round.(Int, LinRange(2, maxmultiple, numthreads+1))
            divide[end] += 1     
            pieces = [divide[k]*i:i:(divide[k+1]-1)*i for k in 1:length(divide)-1]
            @threads for piece in pieces
                sieve[piece] .= false # number 2*i, 3*i, ... are not prime.
            end  
        end 
    end 
    return (1:n)[sieve]
end 
=# 

function removemultiples!(sieve::BitVector, piece::UnitRange{Int64}, i::Int)
    thefirstindex = piece[begin]
    thelastindex = piece[end]
    firstmultiple = max(2, div(thefirstindex, i, RoundUp))
    rangetocheck = i*firstmultiple:i:thelastindex
    view(sieve, rangetocheck) .= false 
    return nothing
end 
#=
for n in [100, 1000, 10_000]
    for i in 1:isqrt(n)
        sieve1 = trues(n)
        sieve2 = trues(n)
        piece = 13:67
        removemultiples!(sieve1, piece, i)
        sieve2[2*i:i:end] .= false 
        @show n i 
        @test sieve1[piece] == sieve2[piece]
    end 
end 
=# 
# start julia with julia --threads=auto or julia --threads=14 # or another number. 
function primesievemultithreaded(n::Integer)
    numthreads = min(nthreads() - 1, n - 2) 
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    # divide work 
    divide = round.(Int, LinRange(2, n, numthreads+1))
    divide[end] += 1     
    pieces = [divide[k]:(divide[k+1]-1) for k in 1:length(divide)-1]        
    for i in 2:isqrt(n)
        if sieve[i] # i is prime 
            for piece in pieces
                removemultiples!(sieve, piece, i)
            end  
        end 
    end 
    return (1:n)[sieve]
end 


function primesieveskippingevens(n::Integer)
    sieve = trues((n-1)÷2) # entry i represents number 2*i+1
    imax = (isqrt(n)-1)÷2
    for i in 1:imax
        if sieve[i] # k=2*i+1 is prime 
            view(sieve, 3*i+1:2*i+1:length(sieve)) .= false # number 3*k, 5*k, ... are not prime (remember that 2*k, 4*k, etc are even and already left out). 
        end 
    end 
    return vcat(2, (3:2:n)[sieve]) # 
end 

# startjulia with julia --threads=auto or julia --threads=14 # or another number. 
function primesieveskippingevensmultithreaded(n::Integer)
    numthreads = nthreads()-1
    sieve =  fill(true, (n-1)÷2) # trues((n-1)÷2) # entry i represents number 2*i+1
    imax = (isqrt(n)-1)÷2
    for i in 1:imax
        if sieve[i] # k=2*i+1 is prime 
            # divide work 
            largest = (n÷(2*i+1)-1)÷2 # 2*largest+1 is the largest odd multiple.
            divide = round.(Int, LinRange(1, largest, numthreads+1))
            divide[end] += 1
            pieces = [(2*divide[index]*i+divide[index]+i):2i+1:(2*(divide[index+1]-1)*i+(divide[index+1]-1)+i) for index in 1:length(divide)-1]
            # testdisjointness(pieces)
            # give each thread a piece of work. 
            @threads for piece in pieces
                view(sieve, piece) .= false
            end 
        end 
    end 
    return vcat(2, (3:2:n)[sieve]) # 
end 


function f(k)
    if isodd(k)
        return 3k+2
    else
        return 3k+1
    end 
end 
function primesieveskippingevensandmultiplesof3(n::Integer)
    # all prime numbers ≥ 5 are of the from 6k-1 or 6k+1. 
    # all numbers ≥ 5 that are 2, 3, 4 modulo 6 = 2*3 are not prime and can be left out.     
    numelements = (n-1)÷6+(n-5)÷6+1 # all numbers > 5 and ≤ n that are 1 or 5 modulo 6.
    sieve = trues(numelements)     
    # when i is odd, sieve[i] represents number 3i+2
    # when i is even, sieve[i] represents number 3i+1. 
    kmax = (isqrt(n)-1)÷3
    for k in 1:kmax
        if sieve[k]  
            if isodd(k) # 3k+2 is prime.
                view(sieve, 5k+3:6k+4:length(sieve)) .= false 
                view(sieve, 7k+4:6k+4:length(sieve)) .= false
            else # k is even, 3k+1 is prime. 
                view(sieve, 5k+1:6k+2:length(sieve)) .= false 
                view(sieve, 7k+2:6k+2:length(sieve)) .= false 
            end  
        end 
    end 
    allprimes = Vector{Int}(undef, 2+sum(sieve)) # allocate memory. 
    allprimes[1] = 2
    allprimes[2] = 3
    vectorindex = 3
    for index in eachindex(sieve)
        if sieve[index]
            allprimes[vectorindex] = f(index)
            vectorindex += 1
        end 
    end 
    return allprimes 
end 

# test whether my function is correct.
# assuming Primes is well implemented. 
#maxnum = 999_999_999
#trueprimes = primes(maxnum)

#=
testnumberset = vcat(5:10, 90:110, 1000:1010, 100_000, 1000_000_000 - 10: 1000_000_000 + 10)
for maxnum in testnumberset # geslaagd. 
    trueprimes = primes(maxnum)
    @test primesieve(maxnum) == trueprimes;
end
println("finished with primesieve") # geslaagd. 
for maxnum in testnumberset # geslaagd. 
    trueprimes = primes(maxnum)
    @test primesievemultithreaded(maxnum) == trueprimes;
end
println("finished with primesievemultithreaded") # geslaagd 
for maxnum in testnumberset # geslaagd. 
    trueprimes = primes(maxnum)
    @test primesieveskippingevens(maxnum) == trueprimes;
end 
println("finished with primesieveskippingevens") # geslaagd.
for maxnum in testnumberset
    @show maxnum
    local trueprimes = primes(maxnum)
    @test primesieveskippingevensmultithreaded(maxnum) == trueprimes;
end 
println("finished with primesieveskippingevensmultithreaded")
for maxnum in testnumberset # geslaagd 
    trueprimes = primes(maxnum)
    @test primesieveskippingevensandmultiplesof3(maxnum) == trueprimes;
end 
println("finished with primesieveskippingevensandmultiplesof3") # geslaagd. 
=#

# calculate number of passes
function countpasses(func)
    maxnum = 1000_000
    count = 0
    start = time()
    while time() - start < 5.0
        func(maxnum)
        count += 1 
    end 
    return count 
end 

function main()
    for func in [primesieve, primesievemultithreaded, primesieveskippingevens, primesieveskippingevensmultithreaded, primesieveskippingevensandmultiplesof3]
        numpasses = countpasses(func)
        println(func, " ", numpasses)
    end 
end 

main()
#=
@test primesieve(maxnum) == trueprimes;
@test primesievemultithreaded(maxnum) == trueprimes;
@test primesieveskippingevens(maxnum) == trueprimes;
@test primesieveskippingevensmultithreaded(maxnum) == trueprimes;
@test primesieveskippingevensandmultiplesof3(maxnum) == trueprimes;


# test the speed. 
maxnum = 1000_000_000
@btime primesieve(maxnum);
@btime primesievemultithreaded(maxnum);
@btime primesieveskippingevens(maxnum);
@btime primesieveskippingevensmultithreaded(maxnum);
@btime primesieveskippingevensandmultiplesof3(maxnum);
=#

# misschien wordt het probleem veroorzaakt dat twee threads dezelfde byte in the bitarray benaderen en dat het dan fout gaat. 