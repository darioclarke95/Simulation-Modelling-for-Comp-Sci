import numpy as np

def upper_triang_form(A,B):
    """
    Performs transformations to convert A to upper triangular form and
    performs the same transfomations on B

    Returns the resulting matrices
    """
    #I don't want to change the original matrices
    A, B = A.copy(), B.copy()
    Arows, Acols = A.shape

    #Get in upper triangular form
    n = Arows

    for i in range(n):
        #What if the value here is 0?
        #need to swap rows
        #pivot = float(A[i,i])
        col = A[i:,i]
        pRow = -1

        for j in range(len(col)):
            if col[j] != 0.0:
                pRow = j
                break

        if pRow == -1:
            pRow = 0

        pRow += i
        #swap the rows
        A[[i,pRow]] = A[[pRow,i]]

        pivot = A[i,i]

        if i < (n-1):
            for j in range(i+1,n):
                val = A[j,i]
                if val != 0:
                    A[j] = ((A[j] / val) * (-pivot)) + A[i]
                    B[j] = ((B[j] / val) * (-pivot)) + B[i]
    return A, B

def diagonalise(A,B):
    """
    Diagonalises A and performs the same transfomations on B

    Returns the resulting matrices
    """

    #I don't want to change the original matrices
    A, B = A.copy(), B.copy()
    #start from last row and woek up
    #last element is the pivot
    #go through column by column and do eliminations
    Arows, Acols = A.shape
    n = Arows

    for i in range(n-1,0,-1):
        pivot = A[i,i]
        for j in range(i-1, -1,-1):
            val = A[j,i]
            if(val != 0):
                A[j] = ((A[j] / val) * (-pivot)) + A[i]
                B[j] = ((B[j] / val) * (-pivot)) + B[i]
    return A,B

def guassian_jordan_elim(A,B):
    """
    Accepts two matrices (A and B), such that (A|B) is the augmented matrix
    representing the system of equations.  A is an n*n matrix and B is a column
    matrix with n rows.
    Example: given the system of equations 2x + 3y = 7, x + y = 1, A = [[2,3],[1,1]]
    and B = [[7],[1]]


    Returns array with the solutions, for non-singular matrices it returns an
    empty list.  likewise for non-square matrices
    """

    #Make sure its a square matrix
    #Make sure its non-singular
    #Diagonalise Augmented Matrix

    #Validate A
    #Is it an nxn matrix?
    Arows, Acols = A.shape
    Brows, Bcols = B.shape
    n= Arows
    if Arows != Acols:
        return []

    #Validate B
    #does it have N rows
    if Brows != Arows:
        return []

    #Make sure A is non-singular
    if np.linalg.det(A) == 0:
        return []

    ATri, BTri = upper_triang_form(A,B)
    ADiag, BDiag = diagonalise(ATri,BTri)
    #Calculate final results
    result = []
    for i in range(n):
        val = BDiag[i]/ADiag[i,i]
        result.append(round(val.item(),1))

    return result

A = np.array([[1.0,1.0,1.0], [-1.0,-1.0,3.0], [-1.0,-1.0,-1.0]])
B = np.array([[0.0],[3.0],[2.0]])
print(guassian_jordan_elim(A,B))
