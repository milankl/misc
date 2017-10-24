## Do prime factorization
import numpy as np
    
def isprime(p):
    """ returns boolean whether p > 1 is prime or not based on fermat test."""
    if p == 2 or p == 3:
        return True
    
    elif (p % 2 == 0) or (p % 3 == 0):
        return False
    
    else:
        # fermat test
        i = 0
        k = 8   # amount of randomly picked numbers at max
        a = 2   # always do the fermat test with base 2
        while i < k:
            if a**(p-1) % p == 1:   # fermat test successful
                i += 1
            else:                   # fermat test failed, hence no prime
                return False
                
            # randomly pick new base
            a = np.random.randint(2,p-1)
                
        return True

def prime_factors(N):
    """ returns the prime factors of a given number N. Divides iteratively by an increasing sequence of prime numbers. """
    i = 2
    factors = []
    
    if isprime(N):  # initial prime number test
        factors.append(N)
        return factors

    while True: 
        while N/i % 1 == 0:
            factors.append(i)
            print(i)
            N = int(N/i)
            if N == 1:
                return factors
                
            if isprime(N):  # checks whether the new N might be prime
                factors.append(N)
                return factors
    
        i += 1
        # to find the next prime number
        while isprime(i) == False: 
            i += 1
           
def primelist(N):
    """ returns a list of the first N prime numbers. """
    
    v = [2] # first prime number
    i = 3   # start with prime 3
    
    # loop through all odd numbers and check whether they are prime based on isprime function
    while len(v) < N:
        if isprime(i):
            v.append(i)
            
        i += 2
    
    return v
    

        