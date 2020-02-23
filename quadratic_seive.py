import math


def gcd(a, b):
    """ Finds the greatest common divisor
    between two numbers using the Euclidean Algorithm.
    :param a: The first number.
    :param b: The second number.
    :return The greatest common divisor between the two numbers.
    """
    # Base case for the recursive call.
    if a == 0:
        return b

    return gcd(b % a, a)


def get_primes_below(n):
    """ Finds the list of primes below a
    certain value n using the Sieve of Eratosthenes.
    :param n: The upper bound for the list of primes.
    :return: A list of primes below n.
    """
    # There are no primes less than 2.
    if n < 2:
        return []

    is_prime = [True] * n  # Initially set everything to be 'prime'.
    is_prime[0] = False    # 0 is not prime.
    is_prime[1] = False    # 1 is not prime.

    # This is the meat of the Sieve of Eratosthenes.
    # Going from i = 2 to ceil(sqrt(n)), if i is prime, we know that all multiples
    # of i (2i, 3i, ...) are not prime, so mark them as false. We do NOT
    # mark i itself. Furthermore, it is sufficient to start marking false
    # at i^2 because all the smaller multiples of i have already been marked
    # at that point.
    for i in range(2, math.ceil(math.sqrt(n))):
        if is_prime[i]:
            for j in range(i**2, n, i):
                is_prime[j] = False

    returning_prime_list = []  # The list of primes that will be returned.

    # We go through the is_prime list and if a given index is still true,
    # we know that number is truly prime. Append the index to the return list.
    for i in range(n):
        if is_prime[i]:
            returning_prime_list.append(i)

    return returning_prime_list
