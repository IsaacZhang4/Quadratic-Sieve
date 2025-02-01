import QuadraticSieve



if __name__ == "__main__":
    # Below I have listed sample values to factorize while printing intermediate outputs to terminal
    ## Values of n should be odd, composite, and not a perfect power
    ### DISCLAIMER: The last sample number takes around 1-2 minutes to test, however, B ,may also be lowered manually to reduce time (CHECK QuadraticSieve.py)
    n_values = [
        3215031751,  # Comment out if testing n5 or individually
        9912409831,  # Comment out if testing n5 or individually
        37038381852397,  # Comment out if testing n5 or individually
        341550071728321,  # Comment out if testing n5 or individually
        #31868712526338419047  # Comment out if testing n1-n4 or individually
    ]

    for n in n_values:
        print(f'Factorization of {n} is {QuadraticSieve.QS(n)}')  # Non-Trivial Factorization Output
