import numpy as np
import math


eps = 1e-12

def read_input(filename):
    f = open(filename, 'r')

    p = int(f.readline())
    n = int(f.readline())

    temp = p
    if p > n:
        temp = p - (p - n)

    a = list()
    for i in range(temp):
        temp = [float(elem) for elem in f.readline().strip('\n').split(' ')]
        for j in range(i+1):
            a.append(temp[j])

    return p, n, a


def read_normal(filename):
    f = open(filename, 'r')

    p = int(f.readline())
    n = int(f.readline())

    a = list()
    for i in range(p):
        temp = [float(elem) for elem in f.readline().strip('\n').split(' ')]
        a.append(temp)

    return p, n, a


# def get_matrix_pos(pos, length):
#     k = 0
#     for i in range(1, length):
#         if (i*(i+1))/2 > pos:
#             k = i - 1
#             break
#     return k, int(pos - (k*(k+1))/2)


def v(i, j):
    return int((i*(i+1))/2 + j)


def get_pq(a, n):
    curr = -99999
    i_m = 0
    j_m = 0
    for i in range(n):
        for j in range(i+1):
            new = abs(a[v(i, j)])
            if i != j and new > curr:
                curr = new
                i_m = i
                j_m = j
    
    return i_m, j_m


def get_tetha(a, p, q):
    alpha = (a[v(p, p)] - a[v(q, q)]) / 2*a[v(p, q)]
    
    t = 0
    if alpha >= 0:
        t = -alpha + math.sqrt(alpha**2 + 1)
    else:
        t = -alpha - math.sqrt(alpha**2 + 1)
    
    c = 1 / math.sqrt(1 + t**2)
    s = t / math.sqrt(1 + t**2)

    return c, s


def is_diagonal(a, n):
    for i in range(n):
        for j in range(i+1):
            if i != j and a[v(i, j)] != 0:
                return False
    return True


def rotation_matrix(n, p, q, c, s):
    r = np.identity(n)
    
    r[p][p] = r[q][q] = c
    r[p][q] = s
    r[q][p] = -s

    return r


def multiply(m, a, n):
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                # daca e deasupra diagonalei principale, apelam la simetrie
                if j > k:
                    res[i][j] += m[i][k] * a[v(j, k)]
                else:
                    res[i][j] += m[i][k] * a[v(k, j)]
    return res

def recalc_a(a, n, p, q, c, s):
    # b = list()
    # b[:] = a[:]

    # for j in range(n):
    #     if j != p and j != q:
    #         b[v(p, j)] = c*a[v(p, j)] + s*a[v(q, j)]
    
    # for j in range(n):
    #     if j != p and j != q:
    #         b[v(q, j)] = -s*a[v(j, p)] + c*a[v(q, j)]
    #         b[v(j, q)] = -s*a[v(j, p)] + c*a[v(q, j)]

    # for j in range(n):
    #     if j != p and j != q:
    #         b[v(j, p)] = a[v(p, j)]

    # b[v(p, p)] = a[v(p, p)] + t*a[v(p, q)]
    # b[v(q, q)] = a[v(q, q)] - t*a[v(p, q)]
    # b[v(p, q)] = 0.0
    # b[v(q, p)] = 0.0

    # for i in range(len(b)):
    #     if -eps < b[i] < eps:
    #         b[i] = 0.0

    # return b

    r = rotation_matrix(n, p, q, c, s)
    ra = np.array(multiply(r, a, n))
    rart = np.dot(ra, r.T)

    res = list()
    for i in range(n):
        for j in range(i+1):
            res.append(rart[i][j])

    return res


def recalc_u(u, n, p, q, c, s):
    r = rotation_matrix(n, p, q, c, s)
    u = np.dot(u, r.T)
    # cu = np.zeros((n,n))
    # cu[:,:] = u[:,:]

    # for i in range(n):
    #     u[i][p] = c*u[i][p] + s*u[i][q]
    
    # for i in range(n):
    #     u[i][q] = -s*u[i][p] + c*u[i][q]

    return u


def jacobi(a, n):
    k = 0
    kmax = 10000
    u = np.identity(n)
    p, q = get_pq(a, n)
    c, s = get_tetha(a, p, q)

    while (not is_diagonal(a, n)) and k <= kmax:
        a = recalc_a(a, n, p, q, c, s)
        u = recalc_u(u, n, p, q, c, s)
        p, q = get_pq(a, n)
        c, s = get_tetha(a, p, q)
        k += 1

    return a, u


def diag(a, n):
    d = list()
    for i in range(n):
        d.append(a[v(i, i)])
    return d


def get_lambdas(values, n):
    res = np.identity(n)
    for i in range(n):
        res[i][i] = values[i]
    return res


def lu_decomp(a):
    n = a.shape[0]
    u = a.copy()
    l = np.eye(n, dtype=np.double)

    for i in range(n):
        factor = u[i+1:, i] / u[i, i]
        l[i+1:, i] = factor
        u[i+1:] -= factor[:, np.newaxis] * u[i]

    return l, u


def cholesky(a):
    a = np.array(a)
    l = np.linalg.cholesky(a)
    print(l, l.T)
    next_a = np.dot(l.T, l)

    while np.linalg.norm(next_a - a) > eps:
        l = np.linalg.cholesky(a)
        next_a = np.dot(l.T, l)
        a = next_a

    return next_a


def get_si(s, p, n):
    si = np.zeros((n, p))
    for i in range(len(s)):
        si[i][i] = 1 / s[i]
    return si


if __name__ == '__main__':
    p, n, a = read_input('input.txt')
    pn, nn, an = read_normal('input.txt')

    if p == n:
        # JACOBI METHOD
        values, vectors = jacobi(a, n)

        values = diag(values, n)
        print('Values:')
        print(values)
        print()

        print('Vectors:')
        for i in range(n):
            print(vectors[:,i])
        print()
        
        au = multiply(vectors.T, a, n)
        uv = get_lambdas(values, n)

        print(au)
        print()
        print(uv)
        print()

        print('Error:', np.linalg.norm(au - uv))

        # CHOLESKY METHOD
        print(an)
        try:
            chol = cholesky(an)
        except:
            print('Something went wrong')
    elif p > n:
        # SINGULAR VALUE DECOMPOSITION
        an = np.array(an)
        print(an)
        print()

        u, s, vh = np.linalg.svd(an)

        print('Singular values:', s)
        print()

        print('Matrix rank:', len(s))
        print()

        print('Condition number:', max(s) / min(s))
        print()

        print('Moore-Penrose pseudo-inverse:')
        si = get_si(s, pn, nn)
        ai = np.dot(np.dot(vh, si), u.T)
        print(ai)
        print()

        aj = np.dot(np.linalg.inv(np.dot(an.T, an)), an.T)
        print(aj)
        print('Norm:', np.linalg.norm(ai - aj, 1))