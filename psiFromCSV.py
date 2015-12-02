import numpy as np

def psiFromCSV(fileName,tidyUp = False):
    """
    Takes the path of a .csv-file as input. The last column is presumed to
    encode the "n"-vector enumerating the occurrance of different haplotypes;
    The other columns are presuemd to encode haplotypes.

    If the below were the contents of myPhi.csv:
    "
    0, 1, 1, 0, 0, 2, 0, 0, 3, 3, 22
    0, 0, 2, 1, 0, 3, 0, 0, 2, 1, 18
    1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 12
    0, 1, 1, 0, 2, 2, 0, 0, 3, 3, 6
    0, 1, 1, 0, 0, 2, 0, 0, 0, 3, 6
    1, 2, 1, 0, 0, 0, 2, 0, 0, 0, 4
    1, 1, 1, 3, 0, 3, 0, 0, 3, 3, 3
    1, 1, 1, 0, 0, 2, 0, 0, 3, 3, 3
    0, 1, 1, 2, 0, 2, 0, 0, 3, 3, 3
    1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 3
    0, 1, 1, 0, 2, 2, 0, 1, 3, 3, 3
    0, 3, 1, 0, 0, 2, 0, 0, 3, 3, 3
    3, 1, 1, 0, 0, 2, 0, 2, 3, 3, 3
    0, 3, 1, 1, 0, 3, 0, 0, 0, 1, 2
    1, 0, 1, 0, 0, 2, 1, 0, 0, 0, 2
    0, 0, 1, 1, 2, 3, 0, 0, 2, 1, 2
    0, 0, 1, 1, 0, 3, 0, 0, 0, 1, 1
    0, 0, 2, 1, 0, 3, 0, 0, 0, 1, 1
    1, 0, 1, 1, 2, 0, 0, 0, 0, 0, 1
    1, 1, 1, 3, 0, 3, 0, 3, 3, 3, 1
    3, 2, 2, 1, 0, 3, 0, 0, 0, 1, 1
    "
    Then psiFromCSV("myPhi.csv") would return S,n whereby:

    >>> S
    array([[0, 1, 1, 0, 0, 2, 0, 0, 3, 3],
       [0, 0, 2, 1, 0, 3, 0, 0, 2, 1],
       [1, 0, 1, 0, 0, 2, 0, 0, 0, 0],
       [0, 1, 1, 0, 2, 2, 0, 0, 3, 3],
       [0, 1, 1, 0, 0, 2, 0, 0, 0, 3],
       [1, 2, 1, 0, 0, 0, 2, 0, 0, 0],
       [1, 1, 1, 3, 0, 3, 0, 0, 3, 3],
       [1, 1, 1, 0, 0, 2, 0, 0, 3, 3],
       [0, 1, 1, 2, 0, 2, 0, 0, 3, 3],
       [1, 0, 1, 1, 2, 2, 0, 0, 0, 0],
       [0, 1, 1, 0, 2, 2, 0, 1, 3, 3],
       [0, 3, 1, 0, 0, 2, 0, 0, 3, 3],
       [3, 1, 1, 0, 0, 2, 0, 2, 3, 3],
       [0, 3, 1, 1, 0, 3, 0, 0, 0, 1],
       [1, 0, 1, 0, 0, 2, 1, 0, 0, 0],
       [0, 0, 1, 1, 2, 3, 0, 0, 2, 1],
       [0, 0, 1, 1, 0, 3, 0, 0, 0, 1],
       [0, 0, 2, 1, 0, 3, 0, 0, 0, 1],
       [1, 0, 1, 1, 2, 0, 0, 0, 0, 0],
       [1, 1, 1, 3, 0, 3, 0, 3, 3, 3],
       [3, 2, 2, 1, 0, 3, 0, 0, 0, 1]])
    >>> n
    array([22, 18, 12,  6,  6,  4,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  1,
        1,  1,  1,  1])
    """
    raw = np.genfromtxt(fileName, delimiter=',')
    # print raw
    if raw.ndim == 1: # If we only have one line, raw will have only 1 axis
        n = np.array(raw[-1], dtype = int, ndmin = 1)
        S = np.array(raw[:-1], dtype = int,ndmin = 2)
    else:
        n = np.array(raw[:,-1], dtype = int)
        S = np.array(raw[:,:-1], dtype = int)

    if tidyUp:
        S = withoutNullColumns(S)
        S,n = removeDuplicateRows(S,n)

    return S,n

def withoutNullColumns(S):
    rows,columns = S.shape
    nonNullColumns = [i for i in range(columns) if max(S[:,i]) > 0]
    columns_new = len(nonNullColumns)

    if columns_new > 0:
        S_new = S[:,nonNullColumns]
        S_new.shape = (rows,columns_new)
    else:
        S_new = np.zeros((rows,1))

    return S_new


def removeDuplicateRows(S,n):

    rows,columns = S.shape
    row_counts = {}

    rows = S.shape[0]

    for row_index in range(rows):
        row = S[row_index]
        rowAsTuple = tuple(row)
        if rowAsTuple in row_counts:
            row_counts[rowAsTuple] += n[row_index]
        else:
            row_counts[rowAsTuple] = n[row_index]

    rowList = row_counts.keys()
    rowCounts = row_counts.values()

    S_new = np.array(np.r_["0,2",rowList],dtype = int)
    n_vec = np.array(rowCounts,dtype = int)

    return S_new,n_vec
