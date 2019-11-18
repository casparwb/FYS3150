include("functions.jl")

function analytical(T::Real, ssave=false)
    """ Calculates the analytical solution to the quantites
    - Z
    - ⟨E⟩, ⟨E^2⟩, σE, Cv
    - ⟨|M|⟩, ⟨|M^2|⟩, σM, χ

    for 2x2 spin lattice at temperature T
    """

    β = 1.0/T

    Es = [-8, 0, 0, 8, 0, -8]   # possible energy values
    Ms = [4, 2, 0, 0, -2, -4]   # possible magnetization values
    degens = [1, 4, 4, 2, 4, 1] # degenercy for each state

    Z = 2*exp(8*β) + 12 + 2*exp(-8*β) # partition function

    """ mean energy, mean squared energy, standard deviation of energy """
    Emean = 1.0/Z*(-16*exp(8*β) + 16*exp(-8*β))
    E2mean = 1.0/Z*(128*exp(8*β) + 128*exp(-8*β))


    Evar = (E2mean - Emean^2)/(4*T^2) # heat capacity

    """ mean magnetization, mean squared magnetization, std. deviation of magnetization """
    Mmean = 0
    M2mean = 16/Z*(2*exp(8*β) + 1)
    Mabs = 16/Z*(exp(8*β) + 1)

    Mvar = (M2mean - Mmean^2)/(4*T) # magnetic susceptibility

    """ per particle """

    Emean /= 4
    Mmean /= 4
    Mabs /= 4

    expecvalues = [Emean, Evar, Mmean, Mvar, Mabs]

    if ssave
        savepath = string(datapath, "analytical_T=$(Float64(T))", ".jld")
        save(savepath, "expecvalues", expecvalues)
    end

    return expecvalues
end
