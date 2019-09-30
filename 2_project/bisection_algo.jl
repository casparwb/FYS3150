
"""
n: order of matrix
c: n x 1 array with the diagonal elements
b: n x 1 array with sub-diagonal elements, b[1] = 0
beta = n x 1 array with the squares of the sub-diagonal elements
m1, m2: eigvals in [lambda_m1, lambda_m1+1,..., lambda_m2]
eps1: tolerance for computed eigvals precision


function bisect(c, b, n)

    beta = [b[i]*b[i] for i=1:n]

    # xmin and xmax:

    xmin = c[n] - abs(b[n])
    xmax = c[n] + abs(b[n])

    for i = 1:n-1
        h = abs(b[i]) + abs(b[i+1])
        if (c[i] + h) > xmax xmax = c[i] + h end
        if (c[i] - h) < xmin xmin = c[i] - h end
    end

    
"""
