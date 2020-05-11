import quadratic_seive as qs
import math
import time
import argparse


def init_parser():
    """ The good ol' argument parser.
    :return: The command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('n', type=int,
                        help='The number you want to factor.')
    parser.add_argument('sb', type=int,
                        help='The smoothness bound used to calculate the factor base.')
    parser.add_argument('si', type=int,
                        help='The upper limit of the sieving interval [0, si]. ')
    parser.add_argument('-v', '--verbose',
                        help='Prints extra information when running QS.', action='store_true')


    return parser.parse_args()


'''
def optimal_sieving_interval_size(n):
    """ The optimal size of the sieving interval as defined in
    https://www.cs.virginia.edu/crab/QFS_Simple.pdf.
    :param n: The number to be factored
    :return: An optimal sieving interval size.
    """
    return round((math.e**(math.sqrt(math.log(n) * math.log(math.log(n)))))**((3 * math.sqrt(2))/4))
'''

if __name__ == "__main__":

    args = init_parser()

    start = time.time()

    n = args.n
    sb = args.sb
    si = args.si

    primes = qs.get_primes_leq(sb)

    fb = qs.make_factor_base(n, primes)
    if fb is None:
        exit(-1)

    non_modded_nums, sieving_nums, factorizations = qs.get_smooth_nums(n, fb, si, primes, args.verbose)
    if non_modded_nums is None:
        exit(-1)

    matrix = qs.build_matrix(factorizations, fb)

    gauss_matrix, pivot, pivot_col_to_row = qs.fast_guass(matrix)

    dependent_rows = qs.find_dependent_rows(gauss_matrix, pivot, pivot_col_to_row)

    solution = qs.solve(dependent_rows, non_modded_nums, sieving_nums, n)
    if solution is None:
        exit(-1)

    if args.verbose:
        # Print out the factor base.
        print('\nFactor Base:')
        print(fb)

        # Print out the smooth numbers and their prime factorizations
        print('\nSmooth Numbers and Prime Factorization:')
        for i in range(len(non_modded_nums)):
            print('{}^2 ≡ {} (mod {}) -> {}'.format(
                int(math.sqrt(non_modded_nums[i])), sieving_nums[i], n, factorizations[i]))

        # Print out the exponent matrix.
        print('\nExponent Matrix (mod 2):')
        for i, row in enumerate(matrix):
            print('Row {}: {} -> {}^2 ≡ {} (mod {})'.format(
                i, row, int(math.sqrt(non_modded_nums[i])), sieving_nums[i], n))

        # Print out the dependent rows, i.e. rows that add to zero modulo 2.
        print('\nDependent Rows (0 indexed):')
        for dependency in dependent_rows:
            print(dependency)

    # Print out the solution.
    print('\nSolution:')
    print('{} = {} * {}'.format(n, solution[0], solution[1]))

    # Print out the total time taken.
    print('\nTime taken:')
    print(time.time() - start, 'seconds.')
