#""" Unit tests """
using Test


function test_RungeKutta()
    """ Test the Runge Kutta implementation with y'=-y """

    T = 10   # maximum time
    n = 1000 # number of grid points

    func(y, t, p) = -y # function to integrate

    t = range(0, stop=T, length=n+1) # time points

    exact_values = exp.(-t)
    exact_values = reshape(exact_values, (length(exact_values),1))

    u0 = 1 # initial value

    approx_values = solve_RK4(func, u0, T, n, nothing)[1] # compute approximated

    @test exact_values ≈ approx_values atol=1e-8 # test that arrays are the same within tolerance 1e-8
end

"""
julia> test_RungeKutta()
Test Passed
"""
