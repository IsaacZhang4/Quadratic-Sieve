from itertools import chain  # Used for matrix manipulation

def printMatrix(matrix):  # helps display matrix
    for row in matrix:
        print(row)

def transpose(matrix):  # Used for transposing a matrix by swapping its rows with columns.
    # This often makes list comprehension easier
    newMatrix = []  # Array to store the transposed matrix
    # Iterate over the column indices (number of columns = len(matrix[0]))
    for i in range(len(matrix[0])):
        newRow = []  # Initialize a new row for the transposed matrix

        # Collect the elements from the current column across all rows
        for row in matrix:
            newRow.append(row[i])  # Append the i-th element of each row to the new row
        newMatrix.append(newRow)  # Add the new row to the transposed matrix

    return (newMatrix)  # Return the transposed matrix

def buildMatrix(smoothNums, factorBase):  # This generates a matrix of exponent vectors (mod 2) based on previously obtained B-smooth numbers.
    # This matrix will be used for the factorization

    factorizationed = []  # Stores the prime factorizations of each smooth number

    def factor(n, factorBase):  # Used to perform trial division to factorize n using the primes in the factor base
        factors = []  # Array to store the prime factors
        if n < 0:
            factors.append(-1)  # Include -1 for negative numbers
        for p in factorBase:
            if p == -1:   # Make sure to skip -1 during the trial division
                pass
            else:
                while n % p == 0:  # Dividing out all of the factors of p
                    factors.append(p)
                    n //= p
        return factors

    matrix = []  # Use this array matrix to store the exponent vectors (mod 2)
    factorBase.insert(0, -1)  # Account for sign by adding -1 to the factor base

    # This processes each smooth number to build exponent vectors
    for n in smoothNums:
        expVector = [0] * (len(factorBase))  # Initialize the exponent vector to all zeros
        nFactors = factor(n, factorBase)  # Factorize the number using the factor base
        factorizationed.append(nFactors)  # Store the factorization result for display

        # Start counting the exponents of each factor(mod 2)
        for i in range(len(factorBase)):
            if factorBase[i] in nFactors:
                expVector[i] = (expVector[i] + nFactors.count(factorBase[i])) % 2  # Add the count of factor occurrences and take mod 2 for parity

        # This checks for a perfect square that is, if all exponents are even then the vector is zero (mod 2)
        if 1 not in expVector:
            return True, n  # If a zero vector is found, n is a perfect square

        matrix.append(expVector)  # Append the exponent vector to the matrix

    # Printing the factorized matrix for part 1.3
    print("Printing x^2-n factorized matrix")  # Intermediate Output
    printMatrix(factorizationed)

    # Return the transposed matrix if no perfect square is found
    return (False, transpose(matrix))

def gaussElim(matrix):  # This uses a simple version of Gaussian elimination to find the reduced row echelon form of a binary matrix.
    # This Gaussian Elimination algorithm identifies the nullspace of the matrix by locating free columns.

    pivotMarks = [False] * len(matrix[0])  # Tracks whether a column contains a pivot element

    # Use Gaussian Elimination to obtain the RREF
    for i in range(len(matrix)):  # Iterate through all rows
        row = matrix[i]  # Current row

        for num in row:  # Search for a pivot element which is identified by 1 in the row
            if num == 1:
                j = row.index(num)  # Column index of the pivot
                pivotMarks[j] = True  # Mark this column as containing a pivot

                # Eliminate 1s in the same column for other rows
                for k in chain(range(0, i), range(i + 1, len(matrix))):  # Iterate over all other rows
                    if matrix[k][j] == 1:  # Check if a 1 is found in the pivot column
                        for i in range(len(matrix[k])):  # Perform row addition (mod 2)
                            matrix[k][i] = (matrix[k][i] + row[i]) % 2
                break  # Move to the next row after processing the pivot column

    # Transpose the matrix to simplify access to columns
    matrix = transpose(matrix)

    solnRows = []  # Stores free rows (corresponding to free columns) and their indices
    # Identify free columns (no pivots)
    for i in range(len(pivotMarks)):
        if pivotMarks[i] == False:  # If the column does not contain a pivot, it's free
            free_row = [matrix[i], i]  # Collect the free row and its index
            solnRows.append(free_row)

    # Return results or handle failure
    if not solnRows:  # If no free rows (solutions) are found
        return ("We could find no solution, must get more smooth numbers.")

    print(f'We found {len(solnRows)} potential solutions') # For debugging
    # Return potential solutions, pivot marks, and modified matrix
    return solnRows, pivotMarks, matrix

def solveRows(solnRows, matrix, pivotMarks, K=0):  # Solves for a row corresponding to a free variable in the e matrix.
    # Determines dependent rows based on the free row and its pivot columns then builds the solution vector.
    solnVec, indices = [], []  # Initialize solution vector and list of pivot column indices
    free_row = solnRows[K][0]  # Extract the K-th free row from solnRows, may be multiple K
    # Identify which columns in the free row contain 1s
    for i in range(len(free_row)):
        if free_row[i] == 1:
            indices.append(i)  # Add column indices where the free row has 1s

    # Identify dependent rows based on pivot columns in the free row
    for r in range(len(matrix)):  # Iterate over all rows in the matrix
        for i in indices:  # Check the columns with 1s in the free row
            if matrix[r][i] == 1 and pivotMarks[r]:  # If a row contains 1 in the same column as pivot
                solnVec.append(r)  # Mark this row as dependent and add to the solution vector
                break

    #  Append the index of the current free row to the solution vector
    solnVec.append(solnRows[K][1])  # Add the index of the free row to the solution vector

    # Return the solution vector and the free row with its index
    return (solnVec), solnRows[K]


