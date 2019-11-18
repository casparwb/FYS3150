using Plots, Random
include("functions.jl")
const func = functions

function montecarlo_new(T::Real, N::Real, L::Real, ordered=false)
    """ Perform N Monte Carlo simulations """

    Es = zeros(Float64, 17) # array for storing all possible values of exp(-ΔE/T)

    """ arrays for storing expectation values """
    eE = zeros(Float64, N)
    eE2 = zeros(Float64, N)
    eM = zeros(Float64, N)
    eM2 = zeros(Float64, N)
    eMabs = zeros(Float64, N)

    ΔE = [-8., -4., 0., 4., 8.] # possible energy difference values
    Es[[1, 5, 9, 13, 17]] = [exp(-E/T) for E in ΔE]

    lattice, E, M = func.init(L, T, ordered)

    @inbounds for iter = 1:N

        """ Metropolis """
        @inbounds for j = 2:L, i = 2:L
            x, y = rand(2:L, 2)

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

    """ Store expectation values """
    eE[iter] = E
    eE2[iter] = E*E
    eM[iter] = M
    eM2[iter] = M*M
    eMabs[iter] = abs(M)


    end # end Monte Carlo

    return eE/N, eE2/N, eM/N, eM2/N, eMabs/N
end

function plot_expecvalues(expecvalues)

    labels = ["E", "E^2", "M", "M^2", "Mabs"]
    cycles = collect(1:length(expecvalues[1]))

    for i = 1:length(expecvalues)
        plot(cycles, expecvalues[i], xlabel="Monte Carlo cycle",
            ylabel=labels[i])
        savefig("figures/expec_"*labels[i]*".png")

    end
end

function main(N::Int64, T::Real, L::Real, ordered)
    """
    Studying how many Monte Carlo cycles are required
    in order to reach the most likely state
    N: number of Monte Carlo cycles
    """


    expecvalues = montecarlo_new(T, N, L, ordered)

    plot_expecvalues(expecvalues)

end
