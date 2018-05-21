from itertools import combinations
from itertools import product
load('mat_to_list.sage')


############################################################
####        Conversion Tools
############################################################


def array_to_tsscpp(A):
    """
    array_to_tsscpp(A)

    Takes a TSSCPP array (see Section 6.2 of "Proofs and
    Confirmations" by David M. Bressoud), and returns the
    encoded TSSCPP in the form of a 2n-by-2n "bird's-eye
    view" matrix

    Inputs
    A: An n-by-n integer matrix encoding a TSSCPP array

    Output
    M: a 2n-by-2n bird's-eye view matrix of the TSSCPP
       encoded in A
    """
    n = A.ncols()
    M = matrix(ZZ, 2*n, 2*n)

    # set the top-left quadrant
    M[0:n, 0:n] = n
    for i in [0..n-1]:
        # fill in blocks belonging to ith shell
        for j in [0..n-1]:
            a = A[i,j]
            k = (a-1)/2
            if a != 0:
                M[j,j] = M[j,j] + 1
                for l in [1..k]:
                    M[j+l,j] = M[j+l,j] + 1
                    M[j,j+l] = M[j,j+l] + 1
        # fill in blocks belonging to previous shells
        for l in [0..n-1]:
            for m in [0..n-1]:
                if l <= i-1 or m <= i-1:
                    M[l,m] = M[l,m] + 1

    # set diagonal
    for l in [0..2*n-1]:
        M[l,2*n-1-l] = n

    # set bottom-right quadrant
    for l in [0..n-1]:
        for m in [0..n-1]:
            M[n+l,n+m] = 2*n - M[n-1-l,n-1-m]

    # set remaining off-diagonal entries
    for l in [0..n-2]:
        for m in [0..n-l-2]:
            c = 0
            for r in [0..n-1]:
                if M[2*n-r-1, n+m] <= l:
                    c = c + 1
                else:
                    break
            M[l, n+m] = 2*n - c
            M[n+m, l] = 2*n - c
            M[n-m-1, 2*n-l-1] = c
            M[2*n-l-1, n-m-1] = c
    return(M)

def tsscpp_to_array(M):
    """
    tsscpp_to_array(M)

    Takes a TSSCPP in the form of a bird's-eye view matrix
    M and returns the TSSCPP array that encodes M. Is the
    inverse of function array_to_tsscpp(A)

    Arguments
    M: A 2n-by-2n bird's-eye view matrix of the TSSCPP

    Outputs
    A: An n-by-n TSSCPP array encoding M
    """

    n = M.nrows()/2
    A = Matrix(ZZ, n, n)

    for i in [0..n-1]:
        partial_shell = Matrix(ZZ, n-i, n-i)
        for l in [0..n-i-1]:
            for m in [0..n-i-1]:
                if M[l+i,m+i] >= 2*n - i:
                    partial_shell[l,m] = 1

        for l in [0..n-i-1]:
            if partial_shell[l,l] == 1:
                counter = 1
                for m in [l+1..n-i-1]:
                    if partial_shell[l,m] == 1:
                        counter = counter + 1
                    else:
                        break
                A[i,i+l] = 2*counter - 1
    return(A)

def plot_tsscpp(M):
    """
    plot_tsscpp(M)

    Takes a TSSCPP in the form of a bird's-eye-view matrix M
    and prints a plot of the TSSCPP

    Inputs:
    M: A 2n-by-2n bird's-eye view matrix of the TSSCPP

    Output:
    (None; prints graphics object)
    """

    PPlist = mat_to_list(M)
    PP = PlanePartition(PPlist)
    return(PP.plot(show_box=True, colors=['ghostwhite', 'lightsteelblue', 'slategrey']))

def array_to_nilp(A):
    """
    array_to_nilp(A)

    Takes a TSSCPP array and returns a plot (as a graphics
    object) of the nest of NILP's associated to the TSSCPP

    Inputs:
    A: An n-by-n integer matrix encoding a TSSCPP array

    Output:
    plot_frame: plot of the nest of NILP's
    """

    n = A.nrows()
    if 2*n-2 <= 12:
        step = 2
    else:
        step = 5

    verts = [[0,0], [n-1,2*n-2], [2*n-2,2*n-2]]
    marks = [1 .. 2*n-2]
    labels = ['${0}$'.format(i) if i % step == 0 else '' for i in [1 .. 2*n-2]]

    plot_frame = polygon(verts, color='black', fill=False, gridlines='major',
                         ticks=[marks, marks], tick_formatter=[labels, labels])
    line_color = 'blue'
    line_weight = 2

    for ii in [1 .. n-1]:
        num_entries = 0
        for kk in [ii .. n-1]:
            if A[ii-1, kk] != 0:
                num_entries += 1
            else:
                break

        x_start = n-ii
        y_start = 2*(n-ii)
        line_points = [(x_start, y_start)]

        if num_entries == 0:
            line_points.append((y_start, y_start))
        else:
            for kk in [1 .. num_entries]:
                entry = A[ii-1, ii-1+kk]
                dist = (entry + 1)/2
                x_upper = 2*x_start - (kk-1) - dist
                y_upper = 2*x_start - (kk-1)
                line_points.append((x_upper, y_upper))
                line_points.append((x_upper, y_upper-1))
                if kk == num_entries and x_upper != y_upper-1:
                    line_points.append((y_upper-1, y_upper-1))

        plot_line = line(line_points, color=line_color, thickness=line_weight, alpha=1.0)
        plot_frame += plot_line
    return(plot_frame)


############################################################
####        Creation Tools
############################################################


def powerset_ordered(S):
    """
    Produces the power set (the collection of all subsets)
    of a set S

    Inputs:
    S: A finite set, as a list of integers

    Output:
    powerset: The power set of S, as a list whose elements
              are the subsets of S, also as lists
    """
    powerset = []

    size = len(S)
    for n_elems in [0..size]:
        if n_elems == 0:
            powerset.append([])
        else:
            combs = list(combinations(S, n_elems))
            n_combs = binomial(size, n_elems)
            for i in [0..n_combs-1]:
                combs[i] = list(combs[i])
            powerset.extend(combs)
    return(powerset)

def make_tsscpp_row(n, i):
    """
    make_tsscpp_row(n,i)

    Intermediate generator used in finding all TSSCPP's
    of order n. This function takes the order of a TSSCPP
    and a row number i and uses an iterator whose elements
    are all possibilities for row i+1 of an n-by-n TSSCPP
    array

    Inputs
    n: Order of the TSSCPP
    i: Row number minus 1 (for purposes of looping)

    Generates
    r: Possible row i+1 for a TSSCPP array of order
       n
    """

    if i == n-1:
        r = [0 for j in [1..n]]
        r[i] = 2*n - 2*i - 1
        yield(r)
        return

    trailing_vals = [-1 + 2*j for j in [1..n-i-1]]
    P = powerset_ordered(trailing_vals)
    P.reverse()
    for subset in P:
        subset.reverse()
        r = [0 for j in [1..n]]
        r[i] = 2*n - 2*i - 1
        m = len(subset)
        if m == 0:
            yield(r)
        else:
            for k in [1..m]:
                r[i+k] = subset[k-1]
            yield(r)

def make_all_rows(n):
    """
    make_all_rows(n)

    Intermediate function that combines all lists of
    possible rows for the TSCPP's of order n

    Inputs
    n: Order of the TSSCPP

    Outputs
    all_rows: A list where each element i (i from 0 to n-1)
              is a list of all possible rows for row i+1 of
              a TSSCPP array of order n
    """

    all_rows = []
    for i in [0 .. n-1]:
        all_rows.append(list(make_tsscpp_row(n, i)))
    return(all_rows)

def make_all_tsscpp_arrays(n):
    """
    make_all_tsscpp_arrays(n):

    Creates a list of all TSSCPP arrays of order n

    Inputs
    n: Order of the TSSCPP array

    Outputs
    all_mats: List of all TSSCPP arrays of order n
    """

    all_mats = []
    R = make_all_rows(n)
    raw_mats = list(product(*R))

    for M in raw_mats:
        mat = []
        for i in [0..n-1]:
            mat.append(M[i])

        i=1
        while 1 <= i and i <= n-1:
            j=i
            while i <= j and j <= n-1:
                if mat[i-1][j] != 0 and mat[i][j] == 0:
                    i = n
                    j = n
                    break
                if mat[i-1][j] > mat[i][j]:
                    i = n
                    j = n
                    break
                if i == n-1 and j == n-1:
                    all_mats.append(Matrix(ZZ, mat))
                j = j+1
            i = i+1
    return(all_mats)
