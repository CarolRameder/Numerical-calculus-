import numpy as np
import math

eps = 1e-12


def read_tridiagonal(filename):
    f = open(filename, 'r')
    n = int(f.readline())
    p = int(f.readline())
    q = int(f.readline())

    sp = f.readline()
    a = [float(f.readline().strip('\n')) for i in range(n)]

    sp = f.readline()
    c = [float(f.readline().strip('\n')) for i in range(n-q)]
    
    sp = f.readline()
    b = [float(f.readline().strip('\n')) for i in range(n-p)]

    return a, b, c


def read_free_terms(filename):
    f = open(filename, 'r')
    n = int(f.readline())
    sp = f.readline()

    a = [float(f.readline().strip('\n')) for i in range(n)]

    return a


def is_positive(a):
    for val in a:
        if val <= eps:
            return False
    return True


def calculate_next(x, a, f):
    a_t, b_t, c_t = a
    norm = 0

    for i in range(len(x)):
        temp = x[i]

        if i == 0:
            x[i] = (f[i] - b_t[i]*x[i+1]) / a_t[i]
        elif i == len(x)-1:
            x[i] = (f[i] - c_t[len(c_t)-1]*x[i-1]) / a_t[i]
        else:
            x[i] = (f[i] - c_t[i-1]*x[i-1] - b_t[i]*x[i+1]) / a_t[i]

        norm += (x[i] - temp)**2

    norm = math.sqrt(norm)
    return x, norm


def multiply(a, x):
    a_t, b_t, c_t = a
    
    sol = list()
    for i in range(len(a_t)):
        if i == 0:
            curr = a_t[0]*x[0] + b_t[0]*x[1]
        elif i == len(a_t)-1:
            curr = c_t[i-1]*x[i-1] + a_t[i]*x[i]
        else:
            curr = c_t[i-1]*x[i-1] + a_t[i]*x[i] + b_t[i]*x[i+1]
        sol.append(curr)

    return sol


def write_to_file(a, filename):
    with open(filename, 'w') as f:
        for list in a:
            f.write('%s\n' % list)


def aprox_solution(a, f):
    a1, b1, c1 = read_tridiagonal(a)
    f1 = read_free_terms(f)

    if is_positive(a1):
        xc = [0 for i in range(len(a1))]
        k = 0
        kmax = 10000
        dx = 1

        while eps <= dx <= 1e+8 and k <= kmax:
            xc, dx = calculate_next(xc, (a1, b1, c1), f1)
            k += 1

        if dx < eps:
            ax = multiply((a1, b1, c1), xc)
            write_to_file(ax, 'ex.txt')
            err = max([abs(ax[i] - f1[i]) for i in range(len(a1))])
            print('Aprox. solution error:', err)
        else:
            print('Divergence...')
    else:
        print('System cannot be solved using Gauss-Seidel...')


if __name__ == '__main__':
    print('a1/f1 - ', end='')
    aprox_solution('a/a1.txt', 'f/f1.txt')
    
    print('a2/f2 - ', end='')
    aprox_solution('a/a2.txt', 'f/f2.txt')
    
    print('a3/f3 - ', end='')
    aprox_solution('a/a3.txt', 'f/f3.txt')
    
    print('a4/f4 - ', end='')
    aprox_solution('a/a4.txt', 'f/f4.txt')
    
    print('a5/f5 - ', end='')
    aprox_solution('a/a5.txt', 'f/f5.txt')
