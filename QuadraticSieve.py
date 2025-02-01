import factorizationHelpers
import linearAlgebraHelpers
import math


def QS(n):  # Quadratic Sieve Algorithm
    root, tolerance = int(math.sqrt(n)), 1

    primality, witnesses = factorizationHelpers.miller_rabin(n)  # Check if n is prime via miller rabin algorithm
    if primality:
        return witnesses

    if isinstance(math.sqrt(n), int):
        return factorizationHelpers.sqrtInt(n)

    print(f'Starting the factorization process of {n}')

    '''
    From testing, I have found that the these values of B tend to work nicer than theory would say
    n1 - B = 373
    n2 - B = 397
    n3 - B = around 1350
    n4 - B = 1447
    n5- B = around 6950
    '''
    B, interval = factorizationHelpers.initializeBounds(n)  # B and interval bounds according to theory

    print(f'Starting to generate a factor base of {B} smoothness')
    factorBase = factorizationHelpers.findBase(n, B)  # generate a factor base of B smoothness
    print(f'These are the primes in factor base: {factorBase}')  # Intermediate Output
    print(f'The total K primes in the factor base is: {len(factorBase)}')  # Intermediate Output


    print(f'Starting to look for {B + tolerance} relations of {B} smoothness')
    smoothNums, xValueList, indices = factorizationHelpers.findBSmooth(factorBase, n, interval, root, tolerance)  # finds relations of B smoothness
    print(f'List of K + 2 values of x: {xValueList}')  # Intermediate Output
    print(f'We found {len(smoothNums)} total numbers of {B} smoothness')  # Debugging

    if len(smoothNums) < len(factorBase):  # If not met, increase B or interval
        return ("Not enough smooth numbers. Increase the sieve interval or size of the factor base.")

    print(f'Starting to build the exponent matrix e')
    isSquareStatus, transpMatrix = linearAlgebraHelpers.buildMatrix(smoothNums, factorBase)  # Builds matrix e

    # Case if completely factorable
    if isSquareStatus == True:
        x = smoothNums.index(transpMatrix)
        factor = int(factorizationHelpers.gcd(xValueList[x] + math.sqrt(transpMatrix), n))
        print("Completely factorable, is a square")
        print("Set if exponents and x values not needed?")  # Intermediate Output
        print(f'Value of x is: {int(xValueList[x])}')  # Intermediate Output
        print(f'Value of y is: {int(math.sqrt(transpMatrix))}')  # Intermediate Output
        print(f'GCD(x-y, n) = {factor}')  # Intermediate Output
        print(f'GCD(x+y, n) = {int(n/factor)}')  # Intermediate Output
        print(f'Finished the factorization process of {n}')
        return factor, int(n / factor)

    print(f'Starting Gaussian Elimination Step')
    solnRows, pivotMarks, matrix = linearAlgebraHelpers.gaussElim(transpMatrix)  # use Gaussian Elimination to solve matrix
    solnVec, solnRowsNew = linearAlgebraHelpers.solveRows(solnRows, matrix, pivotMarks, 0)  # start to find perfect square

    print(f'Starting to solve the congruence of squares')
    factor1, factor2, xValue, yValue = factorizationHelpers.solve(solnVec, smoothNums, xValueList, n)  # solves the congruence of squares

    for J in range(1, len(solnRows)):
        if (factor1 == 1 or factor1 == n):  # Check if solution vector works
            print('Solution vector did not work, attempting new vector')
            solnVec, solnRowsNew = linearAlgebraHelpers.solveRows(solnRows, matrix, pivotMarks, J)
            factor1, factor2, xValue, yValue = factorizationHelpers.solve(solnVec, smoothNums, xValueList, n)
        else:
            print("Found nontrivial factors, printing results")
            print(f'Set of exponents e: {solnRowsNew[0]}\n with length of: {len(solnRowsNew[0])}')  # Intermediate Output
            print(f'x values: {xValueList}\n with length of: {len(xValueList)}')
            print(f'value of x is: {xValue}')  # Intermediate Output
            print(f'value of y is: {yValue}')  # Intermediate Output
            print(f'GCD(x-y, n) = {factor1}')  # Intermediate Output
            print(f'GCD(x+y, n) = {factor2}')  # Intermediate Output
            print(f'Finished the factorization process of {n}')
            return factor1, factor2

    return ("No nontrivial factors found")