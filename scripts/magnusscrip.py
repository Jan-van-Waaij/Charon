from math import ceil
def sieve(limit):
    # Adjust the limit to include the number itself in the sieve.
    limit += 1  
    # Pre-allocate an array for each number, assume all are prime to start.
    # (used to efficiently index into the array to mark non-prime numbers)
    primes = [True] * limit
    # Start eliminating non-prime numbers, starting from 2 up to the square root of the limit.
    for num in range(2, int(limit**0.5 + 1)):
        # Check if the number is a prime number (not divisible by any previous prime number)
        if primes[num]:
            # We now know this number is a prime
            _prime = num
            # How many (whole) times does num go into limit?
            num_of_multiples = ceil(limit / _prime) - 2
            # At every multiple of num, that number is obviously NOT prime
            # (e.g. num = 2, multiples = 4, 6, 8, 10, 12, ...)
            primes[_prime*2:limit:_prime] = [False] * num_of_multiples
    # Set the first two indices to False as 0 and 1 are not prime numbers.
    primes[0] = primes[1] = False
    # Generate prime numbers: For each index, yield the number if its corresponding value is True.
    return list((num for num, is_prime in enumerate(primes) if is_prime))