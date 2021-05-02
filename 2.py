import numpy as np
from scipy.linalg import lu
import math


def get_L(A):
    try:
        for p in range(len(A)):
            A[p][p] = math.sqrt(A[p][p] - sum([A[p][j]**2 for j in range(p)]))

            if p < (len(A)-1):
                for i in range(p + 1, len(A)):
                    A[i][p] = (A[i][p] - sum([A[i][j] * A[p][j] for j in range(p)]))/A[p][p]
        return A
    except Exception as e:
        print('Cannot calculate Cholesky decomposition!', e)


def display(A):
    for i in range(len(A)):
        for j in range(len(A)):
            print(A[i][j], ' ', end='')
        print()
    print()


def display_cholesky(A):
    print('L:')
    for i in range(len(A)):
        for j in range(i+1):
            print(f'{A[i][j]} ', end='')
        for k in range(len(A) - (i+1)):
            print('0.0 ', end='')
        print()

    print()

    print('L^T:')
    for i in range(len(A)):
        for k in range(i):
            print('0.0 ', end='')
        for j in range(i, len(A)):
            print(f'{A[j][i]} ', end='')
        print()


def get_determinant(A):
    det = 1
    for i in range(len(A)):
        det *= A[i][i]
    return det**2


def get_y(A, b):
    y = list()
    for i in range(len(A)):
        y.append((b[i] - sum([A[i][j] * y[j] for j in range(0, i)]))/A[i][i])

    return y


def get_x(A, y):
    x = [0 for i in range(len(A))]

    for i in range(len(A) - 1, -1, -1):
        x[i] = ((y[i] - sum([A[j][i]*x[j] for j in range(len(A) - 1, i, -1)]))/A[i][i])

    return x


def get_initial(A, diag):
    for i in range(len(A)):
        A[i][i] = diag[i]
    for i in range(len(A)):
        for j in range(i+1, len(A)):
            A[j][i] = A[i][j]
    return A


def validate_solution(A, x, b):
    A_init = np.array(A)
    x_chol = np.array(x)
    b_correct = np.array(b)

    return np.linalg.norm(A_init.dot(x_chol) - b_correct)


def lu_decomp(A):
    n = A.shape[0]
    u = A.copy()
    l = np.eye(n, dtype=np.double)

    for i in range(n):
        factor = u[i+1:, i] / u[i, i]
        l[i+1:, i] = factor
        u[i+1:] -= factor[:, np.newaxis] * u[i]

    return l, u


def get_inverse(A):
    A_inverse = list()
    for i in range(len(A)):
         e = [0 for i in range(len(A))]
         e[i] = 1
         x = get_x(A, get_y(A, e))
         A_inverse.append(x)

    return np.array(A_inverse).T


def validate_inv(A_chol, A_bibl):
    return np.linalg.norm(A_chol - A_bibl)


def read_from_file(file):
    with open(file) as f:
        b = [float(n) for n in f.readline().strip().split(' ')]
        A = [[float(n) for n in line.strip().split(' ')] for line in f]
    return A, b


def read_from_input():
    n = int(input('Enter matrix size: '))
    
    print('Enter matrix b:')
    b = list()
    for i in range(n):
        b.append(float(input()))

    print('Enter symmetric matrix A:')
    A = list()
    for i in range(n):
        line = list()
        for j in range(n):
            line.append(float(input()))
        A.append(line)

    return A, b


def generate_random_input():
    n = int(input('Enter matrix size: '))
    A = np.random.rand(n, n)*100
    b = np.random.rand(1, n)*100
    for i in range(len(A)):
        for j in range(i, len(A)):
            A[j][i] = A[i][j]
    return A, b[0]


if __name__ == '__main__':
    method = int(input('Enter method: '))
    if method == 1:
        A, b = read_from_file("input.txt")
    elif method == 2:
        A, b = read_from_input()
    else:
        A, b = generate_random_input()

    print(A)
    print(b)

    diagonal = [A[i][i] for i in range(len(A))]

    if np.linalg.det(A) > 0:
        print('Cholesky decomposition:\n')
        A = get_L(A)
        display_cholesky(A)  # A chol
        print()

        print('Determinant of A:', get_determinant(A))  # A chol

        xchol = get_x(A, get_y(A, b))  # A chol
        print('Cholesky solution of Ax=b equation:', xchol)

        A = get_initial(A, diagonal)
        print('Error of computed result:', validate_solution(A, xchol, b), '\n')  # A init
    
        print('LU decomposition:\n')
        L, U = lu_decomp(np.array(A))  # A init
        print('L:\n', L, '\n')
        print('U:\n', U, '\n')

        xlu = get_x(U.T, get_y(L, b))
        print('LU solution of Ax=b equation:', xlu, '\n')

        A = get_L(A)
        A_inv = get_inverse(A)
        print('A^(-1):\n', A_inv, '\n')
        A = get_initial(A, diagonal)
        print('Error of computed inverse:', validate_inv(A_inv, np.linalg.inv(A)))
    else:
        print('Matrix is not positively defined!')
