# On the cluster, Julia is installed at 
# /home/ctools/julia-1.6.7/bin/julia
# You can install Primes, Test and BenchmarkTools as follows:
#using Pkg
#Pkg.add(["Primes", "Test", "BenchmarkTools"])
#julia> timeremovingnonprimes
#12.42199158668518

#julia> timefilteringoutnumbers
#0.4580650329589844
using BenchmarkTools, Base.Threads, Primes, Test, BenchmarkTools, Base.Iterators, SharedArrays
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

#=
function settofalse(sieve, part, locks)
    for index in part 
        if sieve[index] # still true, skip if already false. 
            gotlocked = trylock(locks[index]) 
            if gotlocked # skip if an other thread is already setting it to false.
                sieve[index] = false
                unlock(locks[index])
            end 
        end 
    end 
end 
function settofalse(sieve, part)
    atomic_xchg!.(view(sieve, part), false)
end 
=#

function settofalse(sieve, part)
    for index in part 
        Base.@atomic sieve[index] = false 
    end 
end 
function primesieveothertry(n::Integer)
    numthreads = nthreads()
    sieve = trues(n) #[Atomic{Bool}(true) for _ in 1:n]
    sieve[1] = false # 1 is not prime. 
    N = n÷(numthreads-1)
    for i in 2:isqrt(n)
        if sieve[i] # i is prime 
            parts = partition(2i:i:n, N)
            Threads.@sync for part in parts 
                Threads.@spawn settofalse(sieve, part)  # number 2*i, 3*i, ... are not prime. 
            end 
        end 
    end 
    return (1:n)[sieve]
end 

maxnum = 1000
trueprimes = primes(maxnum)
@test primesieveothertry(maxnum) == trueprimes

maxnum = 1000_000
@btime primesieve($maxnum);
@btime primesieveothertry($maxnum);

function primesieveotheridea(n::Integer)
    numthreads = nthreads()
    sieve = SharedVector{Bool}(n)
    for index in eachindex(sieve)
        sieve[index] = true 
    end 
    sieve[1] = false # 1 is not prime.
    sqrtn = isqrt(n) 
    @sync for i in 2:sqrtn
        if sieve[i] # i is prime 
            nonprimetoberemoved = i + i 
            while nonprimetoberemoved < sqrtn 
                sieve[nonprimetoberemoved] = false 
                nonprimetoberemoved += i 
            end 
            Threads.@spawn view(sieve, nonprimetoberemoved:i:n) .= false 
        end 
    end 
    return (1:n)[sieve]
end 

maxnum = 1000
trueprimes = primes(maxnum)
@test primesieveotheridea(maxnum) == trueprimes

for _ in 1:1000
    @test primesieveotheridea(maxnum) == trueprimes
end 
maxnum = 1000_000_000 
@btime primesieve($maxnum)
@btime primesieveotheridea($maxnum)



primesieve(1000);
primesieveothertry(1000);
@btime primesieve(10_000_000_000);
@btime primesieveothertry(10_000_000_000);

let 
    vec = trues(1000_000_000)
    t1 = Threads.@spawn vec .= false 
    t2 = Threads.@spawn vec .= false 
    wait(t1)
    wait(t2)
    println(sum(vec))
end 





function primesievetest(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    sqrtn = isqrt(n)
    for i in 2:sqrtn
        if sieve[i] # i is prime 
            view(sieve, 2*i:i:sqrtn) .= false # number 2*i, 3*i, ... are not prime. 
        end 
    end 
    return (1:n)[sieve]
end 

primesieve(1000);
primesievetest(1000);
@btime primesieve(1000_000_000);
@btime primesievetest(1000_000_000);


function removemultiples!(parts, firstindices, index, i)
    firstindexmoduloi = firstindices[index] % i
    if firstindexmoduloi == 0 
        view(parts[index], 1:i:length(parts[index])) .= false 
    else 
        view(parts[index], i-firstindexmoduloi+1:i:length(parts[index])) .= false 
    end 
end  
function primesievesmartidea(n::Integer)
    sqrtn = isqrt(n)
    before = time()
    numthreads = nthreads()
    @assert sqrtn ≤ n÷numthreads
    if n % numthreads == 0 
        sizepart =  n÷numthreads
        parts = [trues(sizepart) for _ in 1:numthreads]
        firstindices = [1 + k*sizepart for k in 0:numthreads-1]
        numparts = numthreads 
    elseif n % (numthreads-1) == 0
        sizepart = n÷(numthreads-1) 
        parts = [trues(sizepart) for _ in 1:numthreads-1]
        firstindices = [1 + k*sizepart for k in 0:numthreads-2]
        numparts = numthreads - 1 
    else 
        sizepart = n÷(numthreads-1) 
        parts = [trues(sizepart) for _ in 1:numthreads-1]
        lastpart = trues(n % (numthreads -1))
        push!(parts, lastpart)
        firstindices = [1 + k*sizepart for k in 0:numthreads-1]
        numparts = numthreads
    end 
    tasksvec = Vector{Task}(undef, numparts)
    timesplittinginparts = time() - before 
    before = time()
    parts[1][1] = false # 1 is not prime. 
    for i in 2:isqrt(n)
        if parts[1][i] # i is prime 
            tasksvec[1] = Threads.@spawn view(parts[1], 2i:i:length(parts[1])) .= false 
            for index in eachindex(parts, firstindices)[2:end]
                tasksvec[index] = Threads.@spawn removemultiples!(parts, firstindices, index, i)
            end 
            wait.(tasksvec)
        end 
    end 
    timeremovingnonprimes = time() - before 
    before = time()
    sieve = vcat(parts...)
    primes = (1:n)[sieve]
    timefilteringoutnumbers = time() - before 
    return timesplittinginparts, timeremovingnonprimes, timefilteringoutnumbers, primes 
end 

maxnum = 1000
timesplittinginparts, timeremovingnonprimes, timefilteringoutnumbers, theprimes  = primesievesmartidea(maxnum);
maxnum = 1000_000_000
timesplittinginparts, timeremovingnonprimes, timefilteringoutnumbers, theprimes  = primesievesmartidea(maxnum)

maxnum = 100_000
trueprimes = primes(maxnum)
@test primesieve(maxnum) == trueprimes
@test primesievesmartidea(maxnum)[4] == trueprimes

maxnum = 1000_000_000
@btime primes(maxnum) seconds=60;
@btime primesievesmartidea(maxnum) seconds=60;

primesieve(1000)
primesieveinbounds(1000)
maxnum = 1000_000_000
@btime primesieve($maxnum) # seconds=600;
@btime primesieveinbounds($maxnum) # seconds=600;



function primesievecounttime(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    timeremovingnonprimes = 0.0 
    timefilteringoutnumbers = 0.0 
    for i in 2:isqrt(n)
        if sieve[i] # i is prime 
            before = time()
            view(sieve, 2*i:i:length(sieve)) .= false # number 2*i, 3*i, ... are not prime. 
            after = time()
            timeremovingnonprimes += after - before 
        end 
    end 
    before = time()
    primes = (1:n)[sieve]
    timefilteringoutnumbers = time() - before 
    return timeremovingnonprimes, timefilteringoutnumbers, primes 
end 

maxnum = 1000
trueprimes = primes(maxnum)
_, _, myprimes = primesievecounttime(maxnum)
@test myprimes == trueprimes

maxnum = 1000_000_000
trueprimes = primes(maxnum)
timeremovingnonprimes, timefilteringoutnumbers, myprimes = primesievecounttime(maxnum)
@test myprimes == trueprimes
timeremovingnonprimes
timefilteringoutnumbers


function primesieve2(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    currentindex = 2 
    sqrtn = isqrt(n)
    while currentindex ≤ sqrtn
        findnextprimetask = Threads.@spawn findnext(sieve, currentindex+1)
        removenonprimestask = Threads.@spawn view(sieve, 2*currentindex:currentindex:length(sieve)) .= false
        wait(removenonprimestask)
        currentindex = fetch(findnextprimetask)
    end 
    return (1:n)[sieve]
end 

function primesieve3(n::Integer)
    sieve = trues(n)
    sieve[1] = false # 1 is not prime. 
    currentindex = 2 
    sqrtn = isqrt(n)
    while currentindex ≤ sqrtn
        view(sieve, 2*currentindex:currentindex:length(sieve)) .= false
        currentindex = findnext(sieve, currentindex+1)
    end 
    return (1:n)[sieve]
end 


maxnum = 100_000
trueprimes = primes(maxnum)
@test primesieve1(maxnum) == trueprimes
@test primesieve2(maxnum) == trueprimes
@test primesieve3(maxnum) == trueprimes

maxnum = 1000_000_000
@btime primesieve1($maxnum) # seconds=600;
@btime primesieve2($maxnum) # seconds=600;
@btime primesieve3($maxnum) # seconds=600;



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

function primesievethreaded2(n::Integer, k::Integer) # is bijna voor elke k even langzaam/snel. getest van 9 tot en met 10*nthreads() = 160. 
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

for m in 1:10*nthreads()
    @btime primesievethreaded2(1000_000_000, $m) seconds=600 
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