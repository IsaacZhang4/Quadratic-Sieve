# Quadratic-Sieve
Overview: Takes a number B for the size of the factor base and applies the program to a number n to return a non-trivial factorization of n.
The program also returns many intermediate values, as these values hold significance related to other congruence or factorization phenomenon.
B is calculated by the program based on the user input n, where n is an odd composite number that isn't a perfect power.

Purpose:
A non-trivial factorization of a number is a factorization that is not 1 and the number itself, this makes the factorization much more difficult to breach ub cryptography. 
This kind of factorization is necessary for many cryptography algorithms, such as RSA, where two non-trivial factors made up of large primes are needed.
The factorization and many of the intermediate values from the program are needed to solve complex modular arithmetic problems and are needed to optimize algorithms such as Pollard's Rho, a crypotography algorithm.
