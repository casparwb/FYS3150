module functions


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
    t = @elapsed @CPUelapsed for i = 2:n
        w = a[i]/b[i-1]

        b[i] -= w*c[i-1]
        d[i] -= w*d[i-1]
    end


    # backward substitution
    t += @elapsed @CPUelapsed v[n] = d[n]/b[n]

    t += @elapsed @CPUelapsed for i = n-1:-1:1
        v[i] = (d[i] - c[i]*v[i+1])/b[i]
    end

    return v, t

end


function exact(x::Float64)
    """ Analytical solution """
    return 1 - (1 - exp(-10))*x - exp(-10*x)
end


function plotting()

    n1 = 10
    n2 = 100
    n3 = 1000

    x1 = range(0, stop=1, length=n1+2)
    x2 = range(0, stop=1, length=n2+2)
    x3 = range(0, stop=1, length=n3+2)

    y1 = matvecmulti_gen(n1)[1]
    y2 = matvecmulti_gen(n2)[1]
    y3 = matvecmulti_gen(n3)[1]

    # include end-points
    y1 = vcat(0, y1, 0)
    y2 = vcat(0, y2, 0)
    y3 = vcat(0, y3, 0)


    # plot each subplot
    pl1 = plot(x1, y1, title="n=10", legend=false)
    pl1 = plot!(x1, map(exact, x1), legend=false)

    pl2 = plot(x2, y2, title="n=100", label="Numerical")
    pl2 = plot!(x2, map(exact, x2), label="Analytical")

    pl3 = plot(x3, y3, title="n=1000", legend=false)
    pl3 = plot!(x3, map(exact, x3), legend=false)

    # plot final figure and save
    plot(pl1, pl2, pl3, layout=(3, 1), xlabel="x", ylabel="u(x)")

    savefig("1bplot.png")

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

    # initialize source term and b
    d = [hh*exp(-10.0*(i*h)) for i = 1:n]
    b = [1.0/((i+1)/i) for i = 1:n]

    # forward substitution
    t = @elapsed @CPUelapsed  for i = 2:n
        d[i] += d[i-1]*b[i-1]
    end

    # backward substitution
    t += @elapsed @CPUelapsed  v[n] = d[n]*b[n]

    t += @elapsed @CPUelapsed  for i = n-1:-1:1
        v[i] = (d[i+1] + v[i+1])*b[i]
    end

    return v, t

end

function CPU_time(n::Int64)


    t_gen = matvecmulti_gen(n)[2]
    t_spec = matvecmulti_spec(n)[2]

    @printf("General algorithm time with %d grid points: %lf ms\n", n, t_gen*1000)
    @printf("Special algorithm time with %d grid points: %lf ms\n", n, t_spec*1000)
    @printf("Relative difference: %lf\n", (t_gen/t_spec))
    @printf("\n")

    return t_gen, t_spec

end


function rel_error(log10n::Int64)
    """log10n: log10 of # grid points"""


    n = 10^log10n
    h = 1.0/(n+1)
    x = range(h, stop=h*n, length=n)

    # calculate the numerical and analytical solutions
    v = matvecmulti_gen(n)[1]
    u = map(exact, x)

    # calculate the absolute difference between each element
    diff = (v .- u) ./ u


    abs_value = broadcast(abs, diff)
    epsilon = maximum(broadcast(log10, abs_value))

    return epsilon
end

function plot_error(plot_, n=1:8)

    err = zeros(length(n)) # array for storing error values

    #calculate, store, and print error values
    for (index, value) in enumerate(n)
        err[index] = rel_error(10^value)
        @printf("n = 10^%d\n", value)
        @printf("log10 relative error: epsilon = %lf\n", err[index])
        @printf("\n")
    end

    # plotting
    if plot_
        x = n
        plot(x, err, xlabel="log10 n", ylabel="log10 max epsilon",
            title="Maximum error as a function of # grid points",
            legend=false)
        scatter!(x, err, legend=false)
        savefig("max_error.png")
    end
end


function compare_LU(n::Int64, julia_tridig=false)
    """A function for timing the built-in lu-factorization function in Julia' linear algebra module

       if julia_tridig is set to true, the matrix A will be constructed using the Tridiagonal-function
       from the linear algebra module, which is more memory efficient."""

    # construct A using LinearAlgebra.Tridiagonal
    if julia_tridig
        dl = [-1 for i = 1:n-1]     # lower and upper diagonal elements

        d = [2 for i = 1:n]         # diagonal elements

        A = Tridiagonal(dl, d, dl)  # create the tridiagonal matrix
    else
    # build the tridiagonal matrix manually
        A = zeros(n, n)
        for i = 1:n
            for j = 1:n
                if (j == i-1) || (j == i+1)
                    A[i, j] = -1
                elseif (j == i)
                    A[i,j] = 2
                end
            end
        end
    end

    t = @elapsed @CPUelapsed lu(A)

    @printf("Julia LU factorization time with %d grid point:  %lf ms\n", n, t*1000)
    @printf("\n")
end

end
