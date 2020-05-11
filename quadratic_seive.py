import math
from decimal import Decimal, getcontext


def gcd(a, b):
    """ Finds the greatest common divisor between two numbers using the Euclidean Algorithm.
    Source: https://en.wikipedia.org/wiki/Euclidean_algorithm
    :param a:
    :param b:
    :return The greatest common divisor between the two numbers a and b.
    """
    # Base case.
    if a == 0:
        return b

    return gcd(b % a, a)  # Recursive call.


def legendre_symbol(a, p):
    """ Calculates the Legendre Symbol.
    Source: https://en.wikipedia.org/wiki/Legendre_symbol
    :param a:
    :param p:
    :return:  1, if the variable a is quadratic residue modulo p and a !≡ 0 (mod p).
             -1, if the variable a is a non-quad residue modulo p.
              0, if a ≡ 0 (mod p).
    """
    return a**((p - 1) // 2) % p


def prime_factorization(n, factor_base):
    """ Calculates the prime factorization of a number against a given factor base.
    :param n: The number we want to find the prime factorization of.
    :param factor_base: A pre-calculated factor base.
    :return: If n factors over the factor base, we return a list containing the prime factorization
             formatted as [[base1, exponent1], [base2, exponent2], ...]. If it doesn't, we return None.
    """
    prime_factors = []
    for prime in factor_base:

        exponent = 0
        while n % prime == 0:
            n //= prime
            exponent += 1

        if exponent != 0:
            prime_factors.append([prime, exponent])

        if prime**2 > n:
            if n > 1:
                prime_factors.append([int(n), 1])

            return prime_factors

    return None


def get_primes_leq(n):
    """ Finds the list of primes at or below a certain value n using the Sieve of Eratosthenes.
    Source: https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
    :param n: The upper bound for the list of primes.
    :return: A list of primes below n.
    """
    # There are no primes less than 2.
    if n < 2:
        return []

    is_prime = [True] * (n + 1)  # Initially set everything to be 'prime'.
    is_prime[0] = False          # 0 is not prime.
    is_prime[1] = False          # 1 is not prime.

    primes = []
    # This is the meat of the Sieve of Eratosthenes.
    # Going from i = 2 to ceil(sqrt(n)), if i is prime, we know that all multiples
    # of i (2i, 3i, ...) are not prime, so mark them as false. We do NOT
    # mark i itself. Furthermore, it is sufficient to start marking false
    # at i^2 because all the smaller multiples of i have already been marked
    # at that point.
    for i in range(2, math.ceil(math.sqrt(n)) + 1):
        if is_prime[i]:
            for j in range(i**2, n + 1, i):
                is_prime[j] = False

    # We go through the is_prime list and if a given index is still true,
    # we know that number is truly prime. Append the index to the return list.
    for i in range(n + 1):
        if is_prime[i]:
            primes.append(i)

    return primes


def make_factor_base(n, primes):
    """ Makes a factor base of primes from the following constraints:
    1. prime <= smoothness_bound.
    2. legendre_symbol(n, prime) = 1.
    :param n: The number to be factored.
    :param primes: The list of primes less than or equal to the smoothness bound.
    :return: A list containing the factor base.
    """

    optimal_size = round((math.e**(math.sqrt(math.log(n) * math.log(math.log(n)))))**(math.sqrt(2)/4))
    factor_base = []

    # Check against the Legendre symbol. If it's equal to 1, put it into the factor base.
    for prime in primes:
        if legendre_symbol(n, prime) == 1:
            factor_base.append(prime)

        if len(factor_base) >= optimal_size:
            break

    if len(factor_base) < optimal_size:
        print('Optimal size of factor base not reached. Please increase the smoothness bound.')
        print('The size of the current factor base is {}. '
              'The optimal size is {}. Do you want to continue QS (Y / N)? '.format(len(factor_base), optimal_size))

        ans = input()
        if ans.lower() == 'y':
            return factor_base
        else:
            return None
    else:
        return factor_base


def get_smooth_nums(n, factor_base, sieving_interval, primes, verbose):
    """ This is the main sieving part of the algorithm. We try to find as many smooth numbers
    as the size of the factor base. This is so that when we build the matrix in the next step,
    it's a square matrix.
    :param n: The number to be factored.
    :param factor_base: The factor base.
    :param sieving_interval: The user-chosen maximum limit M of the sieving interval [0, M].
    :param primes: A list of primes less than or equal to the smoothness bound.
    :param verbose: If the verbose option is set, we print out the progress of the sieve.
    :return:
    """

    floor_sqrt_n = math.floor(math.sqrt(n))
    sieving_list = [(x + floor_sqrt_n)**2 % n for x in range(0, sieving_interval)]

    smooth_seiving_nums = []
    smooth_non_modded_nums = []
    factorizations = []
    for i, seive_num in enumerate(sieving_list):
        non_modded_num = (i + floor_sqrt_n)**2
        prime_factors = prime_factorization(seive_num, primes)

        if prime_factors is None:
            continue

        # A list of bases and exponents of each term in the prime factorization.
        bases = [term[0] for term in prime_factors]
        exponents = [term[1] for term in prime_factors]

        # If all of the bases of each term in the prime factorization are also in
        # the factor base, then it is smooth under the factor base.
        is_smooth = True
        for j in range(len(bases)):
            if bases[j] not in factor_base:
                # We remove the bases not in the factor base from the prime factorization for
                # formatting reasons. There's also a case where the term is base^0, which means
                # that the base is a non-issue when comparing against the factor base, meaning
                # we don't want to return false here. We handle this case next.
                prime_factors.remove([bases[j], exponents[j]])

                if exponents[j] != 0:  # If a base in the prime factorization is not in the factor base AND
                    is_smooth = False  # it has a non-zero exponent, then it's not smooth.

        if is_smooth:
            smooth_seiving_nums.append(seive_num)
            smooth_non_modded_nums.append(non_modded_num)
            factorizations.append(prime_factors)

            if verbose:
                print('Found {} out of {} smooth numbers.'.format(
                    len(smooth_seiving_nums), math.ceil(len(factor_base) + (len(factor_base) * 0.1))))

        if len(smooth_seiving_nums) >= len(factor_base) + (len(factor_base) * 0.1):
            return smooth_non_modded_nums, smooth_seiving_nums, factorizations

    print("Not enough relations. Try increasing sieving interval or smoothness bound.")
    return None, None, None


def build_matrix(factorizations, factor_base):
    """ Builds an exponent matrix from the previously calculated smooth numbers
    and their prime factorizations in GF(2).
    :param factorizations:
    :param factor_base:
    :return:
    """
    matrix = []
    for factorization in factorizations:
        prime_factors = factorization

        exponent_vector = [0] * len(factor_base)

        for j, term in enumerate(prime_factors):
            exponent_vector[factor_base.index(term[0])] = term[1] % 2

        matrix.append(exponent_vector)

    return matrix


def fast_guass(mat):
    """ Used in conjunction with find_dependent_rows to
    find rows in the matrix that sum to the zero vector in modulo 2.
    :param mat:
    :return:
    """

    m_row = len(mat)
    n_col = len(mat[0])

    if m_row < n_col:
        print("More Data needed, Enter more rows.")

    pivot = [False] * m_row
    pivot_found = False
    pivot_col_to_row = {}

    for j in range(n_col):
        pivot_found = False
        # Look for pivot
        for i in range(m_row):
            # Pivot Found at row i and column j
            if mat[i][j] == 1:
                pivot[i] = True
                pivot_col_to_row[j] = i
                pivot_found = True
                break

        if pivot_found:

            for k in range(n_col):

                # Pivot row
                if k == j:
                    continue

                if mat[i][k] == 1:
                    for row_index in range(m_row):
                        mat[row_index][k] = (mat[row_index][j] + mat[row_index][k]) % 2

    return mat, pivot, pivot_col_to_row


def find_dependent_rows(mat, pivot, pivot_col_to_row):
    """ Used in conjunction with fast_gauss to
    find rows in the matrix that sum to the zero vector in modulo 2.
    :param mat:
    :return:
    """

    m_row = len(mat)
    n_col = len(mat[0])

    dependencies = []
    for i in range(m_row):
        # Find Dependent Rows
        if not pivot[i]:
            dep_row = mat[i]
            dependency = [i]
            for j, val in enumerate(dep_row):
                if val == 1:
                    dependency.append(pivot_col_to_row[j])

            dependency = [a for a in dependency]
            dependency.sort()

            dependencies.append(dependency)

    return dependencies


def solve(dependent_rows, non_modded_nums, smooth_nums, n):
    """ Calculates the factorization of n.
    :param dependent_rows: A list of rows in the matrix that sum to zero modulo 2.
    :param non_modded_nums: Numbers from the sieve that are do not have modulo n.
    :param smooth_nums: Numbers that factor over the factor base.
    :param n: The number we want factored.
    :return: The factorization if a non-trivial solution is found. Else, it returns None.
    """
    # Precise calculations up to 1000 digits. Might need to change for extremely large n.
    getcontext().prec = 1000

    for dependency in dependent_rows:
        smooth_solutions = [smooth_nums[row] for row in dependency]
        non_modded_solutions = [math.sqrt(non_modded_nums[row]) for row in dependency]

        smooth_solutions_sq = Decimal(1)
        for num in smooth_solutions:
            smooth_solutions_sq *= Decimal(num)

        non_modded_solutions_sq = Decimal(1)
        for num in non_modded_solutions:
            non_modded_solutions_sq *= Decimal(num)

        factor = gcd(non_modded_solutions_sq - smooth_solutions_sq.sqrt(), Decimal(n))

        if factor != 1 and factor != n:
            return [int(abs(factor)), int(abs(n / factor))]

    print("Cannot find a solution. Try changing the smoothness bound or the sieving interval.")
    return None
