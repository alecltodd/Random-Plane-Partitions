def mat_to_list(M):
    """
    mat_to_list(M)

    Converts a matrix type M into a comma-separated
    list of lists; the elements of the outermost list
    are the rows

    Arguments
    M: a matrix
    """

    m = M.nrows()
    n = M.ncols()
    L = []
    for ii in [0 .. m-1]:
        R = []
        for jj in [0 .. n-1]:
            R.append(M[ii, jj])
        L.append(R)
    return(L)
