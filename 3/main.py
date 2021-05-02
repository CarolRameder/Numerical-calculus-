def read_rare(filename):
    f = open(filename, 'r')
    n = int(f.readline())
    sp = f.readline()

    a = [list() for i in range(n)]
    for line in f:
        line = line.replace(' ', '').strip('\n').split(',')
        val = float(line[0])
        i = int(line[1])
        j = int(line[2])

        new_added = False
        k = 0
        for pair in a[i]:
            if j == pair[1]:
                a[i][k][0] += val
                new_added = True
            k += 1

        if not new_added:
            a[i].append([val, j])

    return a


def read_tridiagonal(filename):
    f = open(filename, 'r')
    n = int(f.readline())
    p = int(f.readline())
    q = int(f.readline())

    sp = f.readline()
    a = [float(f.readline().strip('\n')) for i in range(n)]

    sp = f.readline()
    b = [float(f.readline().strip('\n')) for i in range(n-p)]
    
    sp = f.readline()
    c = [float(f.readline().strip('\n')) for i in range(n-q)]

    return a, b, c


def add(a, b):
    b_a, b_b, b_c = b

    # (i, i)
    for val, i in zip(b_a, range(len(b_a))):
        new_added = False
        k = 0
        for pair in a[i]:
            if i == pair[1]:
                a[i][k][0] += val
                new_added = True
            k += 1

        if not new_added:
            a[i].append([val, i])

    # (i, i+1)
    for val, i in zip(b_b, range(len(b_b))):
        new_added = False
        k = 0
        for pair in a[i]:
            if i+1 == pair[1]:
                a[i][k][0] += val
                new_added = True
            k += 1

        if not new_added:
            a[i].append([val, i+1])

    # (i, i-1)
    for val, i in zip(b_c, range(1, len(b_c)+1)):
        new_added = False
        k = 0
        for pair in a[i]:
            if i-1 == pair[1]:
                a[i][k][0] += val
                new_added = True
            k += 1

        if not new_added:
            a[i].append([val, i-1])

    return a


def multiply(a, b):
    res = [list() for i in range(len(a))]
    b_a, b_b, b_c = b

    for line, i in zip(a, range(len(a))):
        for j in range(len(a)):
            summ = 0
            for pair in line:
                # verificam pe care din diagonale s-ar afla o valoare
                if j == pair[1]:
                    summ += pair[0] * b_a[j]

                elif j+1 == pair[1]:
                    summ += pair[0] * b_c[j]

                elif j-1 == pair[1]:
                    summ += pair[0] * b_b[j-1]

            if summ != 0:
                res[i].append([summ, j])

    return res


# sort by index j
def sort_matrix(a):
    for i in range(len(a)):
        a[i].sort(key=lambda x:x[1])
    return a


def compare_results(a1, a2):
    for list1, list2 in zip(a1, a2):
        if list1 != list2:
            return False
    return True


if __name__ == '__main__':
    a_r = read_rare('a.txt')
    a_t, b_t, c_t = read_tridiagonal('b.txt')

    # added = add(a_r, (a_t, b_t, c_t))
    # true_added = read_rare('aplusb.txt')
    # vf = compare_results(sort_matrix(added), sort_matrix(true_added))
    # print('Added successfully:', vf)

    multiplied = multiply(a_r, (a_t, b_t, c_t))
    true_multiplied = read_rare('aorib.txt')
    vf = compare_results(sort_matrix(multiplied), sort_matrix(true_multiplied))
    print('Multiplied successfully:', vf)
