using Plots
using Printf
using LinearAlgebra
using Statistics


function general_solver(n::Int64)

    """
        A function for solving the general case of a tridigional matrix-vector multiplication
        of size n
    """

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


    # forward substitution, t: CPU time used
    t = @elapsed for i = 2:n
        w = a[i]/b[i-1]

        b[i] -= w*c[i-1]
        d[i] -= w*d[i-1]
    end

    # backward substitution, add the time
    t += @elapsed v[n] = d[n]/b[n]

    t += @elapsed for i = n-1:-1:1
        v[i] = (d[i] - c[i]*v[i+1])/b[i]
    end

    return v, t

end


function exact(x::Float64)
    """ Analytical solution to the differential equation"""
    return 1 - (1 - exp(-10))*x - exp(-10*x)
end


function plotting(n::Int64)
    """A function for plotting the numerical vs. analyical solution
    with varying number of grid points"""

    x1 = range(0, stop=1, length=n+2) # x-values
    x2 = range(0, stop=1, length=10*(n+2))
    y = general_solver(n)[1]        # numerical values from the general algorithm
    y = vcat(0, y, 0)                # add the boundary values (concatenate)

    if n == 100
        plt = plot(x1, y, title="n = $(n)",  xlabel="x", ylabel="v(x)", label="Numerical")
        plt = plot!(x2, map(exact, x2), label="Analytical")
    else
        plt = plot(x1, y, title="n = $(n)", xlabel="x", ylabel="v(x)", legend=false)
        plt = plot!(x2, map(exact, x2), legend=false)
    end

    # save figure
    return plt

end


function special_solver(n::Int64)
    """
        Solver for the general case of the matrix-vector multiplication with
        constant values in the tridiagonal matrix
    """


    h = 1.0/(n + 1)          # step size
    hh = 100*h*h
    d = zeros(Float64, n)    # right side of equation
    b = zeros(Float64, n)    # array for storing 1/b[i]
    v = zeros(Float64, n)    # array for storing result

    # initialize source term and 1/b
    d = [hh*exp(-10.0*(i*h)) for i = 1:n]
    b = [1.0/((i+1)/i) for i = 1:n]

    # forward substitution
    t = @elapsed for i = 2:n
        d[i] += d[i-1]*b[i-1]
    end

    # backward substitution
    t += @elapsed v[n] = d[n]*b[n]

    t += @elapsed for i = n-1:-1:1
        v[i] = (d[i+1] + v[i+1])*b[i]
    end

    return v, t

end

function CPU_time(n::Int64)
    """
        Function for timing the two solvers as a function of n.
        Does an average over 10 calculations
    """
    t_gen = 0
    t_spec = 0
    for i = 1:10
        t_gen += general_solver(n)[2]
        t_spec += special_solver(n)[2]
    end

    t_gen /= 10
    t_spec /= 10

    @printf("General algorithm time with %d grid points: %lf ms\n", n, t_gen*1000)
    @printf("Special algorithm time with %d grid points: %lf ms\n", n, t_spec*1000)
    @printf("Relative difference: %lf\n", (t_gen/t_spec))
    @printf("\n")

    return t_gen, t_spec

end


function rel_error(n::Int64)
    """
        Function for calculating the maximum relative error between the general solver
        and the analytical solution
    """
    h = 1.0/(n+1)                    # step size
    x = range(h, stop=h*n, length=n) # array of x-values to calculate the exact value

    # calculate the numerical and analytical solutions
    v = general_solver(n)[1]
    u = map(exact, x)

    # calculate the difference between each element
    diff = (v .- u) ./ u


    # compute the absolute value of the differences and extract the maximum value
    abs_value = broadcast(abs, diff)
    epsilon = maximum(broadcast(log10, abs_value))

    return epsilon
end

function plot_error(n=1:8, savefigs=false)


    err = zeros(length(n)) # array for storing error values

    #calculate, store, and print error values
    for (index, value) in enumerate(n)
        err[index] = rel_error(10^value)
        @printf("n = 10^%d\n", value)
        @printf("log10 relative error: epsilon = %lf\n", err[index])
        @printf("\n")
    end

    # plotting

    x = n
    plot(x, err, xlabel="log10 n", ylabel="log10 max epsilon",
        title="Maximum error as a function of # grid points",
        legend=false)
    scatter!(x, err, legend=false)
    if savefigs savefig("max_error.png") end

end


function compare_LU(n::Int64)
    """
    A function for timing the built-in lu-factorization function in Julia' linear algebra module
    """

    # build the tridiagonal matrix
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

    # average time over 10 lu calls
    t = 0
    for i = 1:10
        t += @elapsed lu(A)
    end

    t /= 10

    @printf("Julia LU factorization time with %d grid points:  %lf ms\n", n, t*1000)
    @printf("\n")

end

function plot_difference()
    """
       A function for visualizing the difference in the approximated
       values from the general and special algorithm
    """

    n = 1:6
    avg_diff = zeros(Float64, length(n))
    for i in n
        v_gen = general_solver(10^i)[1]
        v_spec = special_solver(10^i)[1]
        avg_diff[i] = sum(v_gen .- v_spec)
        avg_diff[i] /= length(avg_diff[i])
    end
    plot(n, avg_diff)
    savefig("difference.png")
end

function main(n::Int64, savefigs=false)


    # make a subplot with n = 10, 100 and 1000
    pl1 = plotting(10)
    pl2 = plotting(100)
    pl3 = plotting(1000)
    plot(pl1, pl2, pl3, layout=(3, 1))
    if savefigs savefig("compare_general.png") end


    CPU_time(n)

    # only perform the LU factorization for n <= 1000
    if n <= 10^3 compare_LU(n) end


end
