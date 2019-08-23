using Plots

function matvecmulti(n)

    h = 1/(n + 1)          # step size
    a = zeros(Float64, n+1)
    b = zeros(Float64, n+1)
    c = zeros(Float64, n+1)

    v = zeros(Float64, n+1)
    b_tild = zeros(Float64, n+1) # right side of equation

    # initialize
    for i = 1:n+1
        a[i] = -1
        b[i] = 2
        c[i] = -1
        b_tild[i] = 100*exp(-10*i)
    end

    b_tild .*= h^2

    a[1] = 0
    c[n] = 0

    c[1] = c[1]/b[1]
    b_tild[1] = b_tild[1]/b[1]

    for i = 2:n
        d = 1.0/(b[i] - c[i-1]*a[i])

        c[i] *= d
        b_tild[i] = (b_tild[i] - a[i]*b_tild[i-1])*d
    end

    v[n] = b_tild[n]

    for i = n-1:-1:1
        v[i] = b_tild[i] - (c[i]*v[i+1])/b[i]

    end

    return v

end

y = matvecmulti(100)

println(y)
