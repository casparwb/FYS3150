using Plots, LinearAlgebra, JLD, Random, StaticArrays, Formatting, Test, Revise
using PlotThemes; theme(:dark)
const figpath = "figures/"
const datapath = "data/"

include("functions.jl")


function montecarlo(N, T, L, ordered=false)

    """
    Perform N Monte Carlo simulations with Metropolis sampling
    Input
    -------
    N:       number of cycles
    T:       temperature
    L:       number of spins in each dimension
    ordered: true/false - whether to start in an ordered state (all 1s)

    Returns:
    ---------
    -  something
    """

    avgs = zeros(Float64, N, 5)

    Es = zeros(Float64, 17)
    ΔEs = [8., 4., 0., -4., -8.]

    Es[[1, 5, 9, 13, 17]] = [exp(-E/T) for E in ΔEs]
    Es = SVector{17}(Es) # static array for speed

    lattice, E, M = init(L, T, ordered) # initial

    """ Perform N Monte Carlo simulations """
    @inbounds for i = 1:N

        """ Metropolis """
        @inbounds for n = 1:L^2
            x, y = rand(2:L+1, 2)

            ΔE = 2*lattice[y, x]*(lattice[y, x+1] +
                                  lattice[y, x-1] +
                                  lattice[y-1, x] +
                                  lattice[y+1, x])

            w = Es[ΔE + 9]
            if rand() <= w

                lattice[y, x] *= -1
                M += 2*lattice[y, x]
                E += ΔE

            end # end if

        end # end metropolis

    """ Update expectation values """
    avgs[i,:] = [E, E*E, M, M*M, abs(M)]


    end # end Monte Carlo
    return avgs
end




function main(N::Real)
    """ Main program """


    """ When is most likely state reached? """
    L = 20
    T1 = 1.0
    T2 = 2.4


    avgs = montecarlo(N, T2, L, false)

    expecvalues = zeros(Float64, 5)

    expecvalues[:] = compute_quants(sum(avgs, dims=1), [T2, i, L])

    plot(1:N, expecvalues[1:N,1], legend=false)

    # E = zeros(Float64, N÷10)
    # Cv = zeros(Float64, N÷10)
    # M = zeros(Float64, N÷10)
    # χ = zeros(Float64, N÷10)
    # Mabs = zeros(Float64, N÷10

end
main(100)
