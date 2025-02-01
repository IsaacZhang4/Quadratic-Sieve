import QuadraticSieve
import math

def gcd(a, b):  # Standard GCD calculation via Euclid's algorithm
    if b == 0:
        return a
    elif a >= b:
        return gcd(b, a % b)
    else:
        return gcd(b, a)


def LJ(m, n):  # Algorithm for Legendre/Jacobi Symbol (m/n)
    # Note, may need to check if n is odd prime or greater than 2,
    # however, this seemed to bug out my code for the first value of n

    m = m % n  # LJ 1.1 - Compute m ≡ m (mod n), so that the new m satisfies 0 ≤ m < n. Put t = 1.
    t = 1

    while m != 0:  # LJ 1.2 - While m ̸ = 0 {
        while m % 2 == 0:  # LJ 1.2.1. While m is even {put m = m/2 and, if n ≡ 3 or 5 (mod 8), then put t = −t.}
            m //= 2
            if n % 8 == 3 or n % 8 == 5:
                t = -t

        m, n = n, m  # LJ 1.2.2. Interchange m and n.

        if m % 4 == 3 and n % 4 == 3:  # LJ 1.2.3. If m ≡ n ≡ 3 (mod 4), then put t = −t.
            t = -t

        m = m % n  # 1.2.4. Compute m ≡ m (mod n), so that the new m satisfies 0 ≤ m < n.}

    # LJ 2.1. If n = 1, then return t
    # LJ 2.2. Else return 0.
    return t if n == 1 else 0

def QC357(a, p):
    # Note this algorithm only works for a prime p ≡ 3, 5, 7 (mod 8)
    # Check if LJ symbol is 1
    if LJ(a, p) != 1:
        return None

    # Case 1: p ≡ 3 or 7 (mod 8)
    if p % 8 == 3 or p % 8 == 7:
        return fast_mod_exp(a, (p + 1) // 4, p)
    # Case 2: p ≡ 5 (mod 8)
    elif p % 8 == 5:
        x = fast_mod_exp(a, (p + 3) // 8, p)
        # Case 2.1: x^2 ≡ a (mod p)
        if fast_mod_exp(x, 2, p) == a % p:
            return x
        # Case 2.2: not congruent
        else:
            return (x * fast_mod_exp(2, (p - 1) // 4, p)) % p
    # Case 3: this algorithm does not apply
    else:
        return tonelli_shanks(a, p)


def fast_mod_exp(a, expo, m):
    res = 1
    # Make usage of bit shift
    while expo > 1:
        if expo & 1:
            res = (res * a) % m
        a = a ** 2 % m
        expo >>= 1
    return (a * res) % m

def sqrtInt(n):  # Use Newton's method to compute the integer square root of a given number n.
    x = n  # Start with an initial guess for the square root (x = n)
    y = (x + 1) // 2   # Compute the first next approximation using Newton's method's formula
    while y < x:  # Continue iterating until the approximations stabilize (y >= x)
        x = y
        y = (x + n // x) // 2  # This progressively refines the estimate with Newton's method's formula
    return x

# CREDIT FOR THIS CODE - bits taken from ZeroBone GitHub and Utkarsh Trivedi from Geeks for Geeks
def tonelli_shanks(a, p):
    # Check that 'a' is a quadratic residue
    if LJ(a, p) != 1:
        return None

    # Special case where p ≡ 3 mod 4
    if p % 4 == 3:
        x = fast_mod_exp(a, (p + 1) // 4, p)
        return x

    # Get q and s st p - 1 = q * 2 ^ s
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1

    # Get z, a quadratic non-residue (QNR)
    z = 2
    while LJ(z, p) != -1:
        z += 1

    m = s
    c = fast_mod_exp(z, q, p)
    t = fast_mod_exp(a, q, p)
    r = fast_mod_exp(a, (q + 1) // 2, p)

    # Keep looping until soln is found
    while t != 0 and t != 1:
        t2i = t
        i = 0
        for i in range(1, m):
            t2i = fast_mod_exp(t2i, 2, p)
            if t2i == 1:
                break
        b = fast_mod_exp(c, 2 ** (m - i - 1), p)
        r = (r * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i

    return r if t == 1 else None

def is_quadratic_residue(a, p):  # Checks if a is a quadratic residue modulo p using the Legendre symbol
    return LJ(a, p) == 1


def find_b(p):  # Finds an integer b such that (b/p) = -1
    b = 2
    while LJ(b, p) != -1:
        b += 1
    return b


def factor_out_two(p_minus_1):  # Factors p - 1 as 2^s * u where u is odd
    s = 0
    u = p_minus_1
    while u % 2 == 0:
        u //= 2
        s += 1
    return s, u


def compute_m(d, f, s, p):  # Computes m such that (d * f^m) ≡ 1 (mod p) using the method described
    m = 0
    for i in range(s):
        g = (d * pow(f, m, p)) ** (2 ** (s - 1 - i)) % p
        if g == p - 1:
            m += 2 ** i
    return m


def QC1(a, p):
    # Solves x^2 ≡ a (mod p) for a prime p ≡ 1 (mod 8)
    # Step 1: Find b such that (b/p) = -1
    b = find_b(p)

    # Step 2: Factor out powers of 2 in p - 1, so p - 1 = 2^s * u with u odd
    s, u = factor_out_two(p - 1)

    # Step 3: Compute d ≡ a^u (mod p) and f ≡ b^u (mod p)
    d = pow(a, u, p)
    f = pow(b, u, p)

    # Step 4: Compute m such that (d * f^m) ≡ 1 (mod p)
    m = compute_m(d, f, s, p)

    # Step 5: Compute the solution x ≡ a^((u+1)//2) * f^(m//2) (mod p)
    x = (pow(a, (u + 1) // 2, p) * pow(f, m // 2, p)) % p
    return x

def miller_rabin(n):
    # Part 0 -> Check if n is odd, stop if not (witness 2)
    if n % 2 == 0:
        return False, f"{n} is composite because it is even (witness: 2)."

    # Part 1 -> Check n for small factors not exceeding log(n)
    limit = int(math.log(n)) + 1
    for i in range(2, limit):
        if n % i == 0:
            return False, f"{n} is composite because it has a small factor {i} (witness: {i})"

    # Part 2 -> Check if n is a prime power
    max_k = int(math.log(n) / math.log(2))
    for k in range(2, max_k + 1):
        root = int(n ** (1 / k))
        if root ** k == n:
            return False, f"{n} is composite because it is a power of {root} (witness: {root})"

    # Part 3 -> Express n - 1 as (2^u) * v by taking out power of 2 in n - 1
    v = n - 1
    u = 0
    while v % 2 == 0:
        v //= 2
        u += 1

    # Part 4 + 5 -> Check for witnesses
    witness_limit = min(2 * (int(math.log(n)) ** 2), n - 2)
    for a in range(2, witness_limit + 1):
        x = fast_mod_exp(a, v, n)  # a^v % n
        if x == 1 or x == n - 1:
            continue

        # Check the values a^(2^r * v) such that 0 <= r <= u-1
        composite = True
        for _ in range(u - 1):
            x = fast_mod_exp(x, 2, n)
            if x == n - 1:
                composite = False
                break
        if composite:
            return False, f"{n} is composite (witness {a})"

    # Part 6 -> No witness found thus n is prime
    return True, f"{n} is a prime"


def generatePrime(n):  # Use the Sieve of Eratosthenes algorithm to generate a list of prime numbers up to n
    if n < 2: # If n is less than 2, there are no primes; return an empty list
        return []

    # Initialize lists to store numbers and their primality status
    nums = []  # List of numbers from 0 to n
    isPrime = []  # Boolean list to indicate if a number is prime

    for i in range(0, n + 1):
        nums.append(i)  # Fill nums with numbers from 0 to n
        isPrime.append(True)  # Assume all numbers are prime initially

    # 0 and 1 are not prime numbers
    isPrime[0] = False
    isPrime[1] = False

    # Iterate through numbers up to half of n (enough to mark all composites)
    for j in range(2, int(n / 2)):  # iterate all the sensible size gaps
        if isPrime[j] == True:
            # Mark multiples of j as non-prime starting from 2*j
            for i in range(2 * j, n + 1, j):
                isPrime[i] = False
    # Collect all numbers marked as prime into the primes list
    primes = []
    for i in range(0, n + 1):  # Adds leftovers
        if isPrime[i] == True:
            primes.append(nums[i])  # Add prime numbers to the list

    # Return the list of primes
    return primes

def quadResidue(a, n):  # This function Determines whether a is a quadratic residue (modulo n) via LJ
    # Note, may need to check if n is odd prime or greater than 2,
    # however, this seemed to bug out my code for the first value of n

    # Legendre symbol determines quadratic residue
    result = LJ(a, n)

    # Interpret the Legendre symbol result
    if result == 1:
        return True  # Case where ` is a quadratic residue (modulo n)
    elif result == -1:
        return False  # Case where a is not a quadratic residue (modulo n)
    else:  # Essentially result == 0
        return None  # Case where a is 0 (mod n) (neither residue nor non-residue)


def initializeBounds(n):  # Based on QS initialization for B, however, the given B seems too small

    Ln = math.exp(math.sqrt(math.log(n) * math.log(math.log(n))))  # L(n) = exp(sqrt((logn*log(logn))
    B = math.ceil(math.pow(Ln, 1/math.sqrt(2)))  # L(n)^(1/sqrt(2)), ceiling operation
    interval = B ** 2
    return int(B), int(interval)


def findBase(n, B):
    # Generates a B-smooth factor base such that n is a quadratic residue (modulo p)
    # Note that the B-smooth factor base has prime numbers <= B such that n is a quadratic residue.

    factorBase = []  # Initialize empty array for storing the primes in factor base.
    primes = generatePrime(B)  # Generates us all the prime numbers that are less than or equal to B

    for p in primes:  # Iterate through the generated primes
        if quadResidue(n, p) == 1:  # Checks if n is a quadratic residue (modulo p)
            factorBase.append(p)  # If n is a quadratic residue, add the prime p to factor base
    return factorBase  # Outputs a list of primes that form the B-smooth factor base


def findBSmooth(factorBase, n, interval, root, tolerance):  # Used for finding B-smooth numbers in a generated sequence via sieving process
    # Note that a B-smooth number is an integer that factors completely into the prime numbers in the given factor base

    def sieve_prep(n, sieve_int, root):  # This is used for generating a sequence from Y(x) = x^2 - N which will start at x = root
        sieve_seq = [x ** 2 - n for x in range(root, root + sieve_int)]
        return sieve_seq

    def get_residues_for_p(p, n):  # This is used for finding residues x such that x^2 ≡ N (mod p), note, uses QC357 or QC1 depending on p mod 8
        if p % 8 == 1:
            return [QC1(n, p)]  # QC1 is used for the case when p ≡ 1 (mod 8)
        elif p % 8 in {3, 5, 7}:
            result = QC357(n, p)  # QC357 is used here to handle other congruence classes
            return [result] if result is not None else []  # Ensure that the result is iterable
        else:
            return []  # No valid residues for any primes that are outside our defined congruence classes

    # Generate the initial sequence for sieving
    sieve_seq = sieve_prep(n, interval, root)  # Sequence: Y(x) = x^2 - n
    sieve_list = sieve_seq.copy()  # Copy the sieve so it is saved for later modification during sieving

    # Handle the prime factor 2 separately (if it is in the factor base)
    if factorBase[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:  # Finds the first even term in the sequence
            i += 1
        for j in range(i, len(sieve_list), 2):  # Process every other term (all even numbers)
            while sieve_list[j] % 2 == 0:  # Account for all powers of 2 from the sequence by removing them
                sieve_list[j] //= 2

    # Process every prime in the factor base (excluding 2)
    for p in factorBase[1:]:
        residues = get_residues_for_p(p, n)  # Find the solutions to x^2 ≡ n (mod p)

        for r in residues:  # Start position based on the residue and iterate through every pth term
            for i in range((r - root) % p, len(sieve_list), p):
                while sieve_list[i] % p == 0:  # This removes factors of p, which also includes powers of p
                    sieve_list[i] //= p

    # Start the identification of smooth numbers and their corresponding x values
    xlist = []  # List to store x-values such that x^2 - n is smooth
    smoothNums = []  # List to store the numbers of B smoothness
    indices = []  # List to store the corresponding indices of smooth numbers in the sieve sequence


    # Start the process of collecting all the B smooth numbers
    for i in range(len(sieve_list)):
        if len(smoothNums) >= len(factorBase) + tolerance:
            break  # Stop early if enough smooth numbers are found to ensure success probability of 2^-tolerance
        if sieve_list[i] == 1:  # Check if the sequence term is B-smooth, that is, it has been fully reduced
            smoothNums.append(sieve_seq[i])  # Add the original sequence value of B-smoothness
            xlist.append(i + root)  # Add the corresponding x-value where x = root + index
            indices.append(i)  # Add the index where the smooth number was found

    # Output the B-smooth numbers, corresponding x-values, and corresponding index
    return smoothNums, xlist, indices

def solve(solnVec, smoothNums, xlist, n): # Solves for factors of n using the solution vector and smooth numbers obtained from the sieve.
    # The function computes the value of y and x, then uses the GCD to find factors of n.
    # Extrac the smooth numbers and their corresponding x values based on the solution vector indices
    solnNums = [smoothNums[i] for i in solnVec]  # Smooth numbers that contribute to the solution
    xNums = [xlist[i] for i in solnVec]  # Corresponding x-values for those smooth numbers

    # Compute the product of the smooth numbers to obtain Asquare
    smoothProduct = 1
    for num in solnNums:
        smoothProduct *= num  # Multiply all smooth numbers to get Asquare

    # Compute the product of x-values to obtain x
    x = 1

    for num in xNums:
        x *= num  # Multiply all x-values to get x

    y = sqrtInt(smoothProduct)  # Take the integer square root of the product of all the smooth numbers to get y

    # Use GCD to find a factor of n and its corresponding factor
    factor1 = gcd(x - y, n)  # GCD of (b - a) and n gives one factor
    factor2 = gcd(x + y, n)  # GCD of (b + a) and n gives the other factor
    return factor1, factor2, x, y  # Return the two factors and x and y