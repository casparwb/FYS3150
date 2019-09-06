
using Printf
#using Plots

function matvecmulti_gen(n::Int64)

    """Generalized matrix-vector multiplication with a tridiagonal matrix"""

    h = 1.0/(n + 1)          # step size
    hh = 100*h*h
    a = zeros(Float64, n)
    b = zeros(Float64, n)
    c = zeros(Float64, n)
    d = zeros(Float64, n)

    v = zeros(Float64, n)

    # initialize

    for i = 1:n
        a[i] = -1.0
        b[i] = 2.0
        c[i] = -1.0
        d[i] = hh*exp(-10.0*(i*h))  # source term
    end

    a[1] = 0.0
    c[n] = 0.0


    # forward substitution
    for i = 2:n
        w = a[i]/b[i-1]

        b[i] -= w*c[i-1]
        d[i] -= w*d[i-1]
    end


    # backward substitution
    v[n] = d[n]/b[n]

    for i = n-1:-1:1
        v[i] = (d[i] - c[i]*v[i+1])/b[i]
    end

    return v

end


function exact(x::Float64)
    """ Analytical solution """
    return 1 - (1 - exp(-10))*x - exp(-10*x)
end

function matvecmulti_spec(n::Int64)
    """
        Matrix-vector multiplication with a tridiagonal Toeplitz matrix
        n: number of grid points
    """


    h = 1.0/(n + 1)          # step size
    hh = 100*h*h
    d = zeros(Float64, n)    # right side of equation
    b = zeros(Float64, n)    # array for storing 1/b[i]
    v = zeros(Float64, n)    # array for storing result

    # initialize source term and 1.0/b
    d = [hh*exp(-10.0*(i*h)) for i = 1:n]
    b = [1.0/((i+1)/i) for i = 1:n]

    # forward substitution
    @timev for i = 2:n
        d[i] += d[i-1]*b[i-1]
    end

    # backward substitution
    v[n] = d[n]*b[n]

    for i = n-1:-1:1
        v[i] = (d[i+1] + v[i+1])*b[i]
    end

    return v

end

d1 = matvecmulti_spec(10000)
d2 = matvecmulti_gen(10000)

for i = 1:100 @printf("%lf | %lf\n", d1[i], d2[i]) end

function rel_error(n::Int64)

    h = 1.0/(n+1)
    x = range(h, stop=h*n, length=n)

    # calculate the numerical and analytical solutions
    v = matvecmulti_gen(n)
    u = map(exact, x)

    # calculate the absolute difference between each element
    diff = (v .- u) ./ u


    abs_value = broadcast(abs, diff)
    epsilon = maximum(broadcast(log10, abs_value))

    return epsilon
end




# n1 = 10
# n2 = 100
# n3 = 1000
#
# x1 = range(0, stop=1, length=n1)
# x2 = range(0, stop=1, length=n2)
# x3 = range(0, stop=1, length=n3)
#
# y1 = matvecmulti_gen(n1)
# y2 = matvecmulti_gen(n2)
# y3 = matvecmulti_gen(n3)
#
# pl1 = plot(x1, y1, title="n=10", legend=false)
# pl1 = plot!(x1, map(u, x1), legend=false)
#
# pl2 = plot(x2, y2, title="n=100", label="Numerical")
# pl2 = plot!(x2, map(u, x2), label="Analytical")
#
# pl3 = plot(x3, y3, title="n=1000", legend=false)
# pl3 = plot!(x3, map(u, x3), legend=false)
#
# plot(pl1, pl2, pl3, layout=(3, 1), xlabel="x", ylabel="f(x)")
# display(plot)
#savefig("1bplot.png")
