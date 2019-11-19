using Random, StaticArrays

function analytical(T::Real, ssave=false)
    """ Calculates the analytical expectation values
        for the mean energy, magnetization, and the heat capacity
        and magnetic susceptibility for a 2x2 spin lattice at temperature T

    Input:
    ----------------
    - T: temperature [kT/J]

    Output:
    ------------------
    - array with [⟨E⟩, Cv, ⟨M⟩, χ, ⟨|M|⟩] per particle (spin)

    """


    β = 1.0/T

    Es = [-8, 0, 0, 8, 0, -8]   # possible energy values
    Ms = [4, 2, 0, 0, -2, -4]   # possible magnetization values
    degens = [1, 4, 4, 2, 4, 1] # degenercy for each state

    Z = 2*exp(8*β) + 12 + 2*exp(-8*β) # partition function

    """ mean energy, mean squared energy, heat capacity """
    Emean = 1.0/Z*(-16*exp(8*β) + 16*exp(-8*β))
    E2mean = 1.0/Z*(128*exp(8*β) + 128*exp(-8*β))


    Cv = (E2mean - Emean^2)/(4*T^2) # heat capacity

    """ mean magnetization, mean squared magnetization, magnetic susceptibility """
    Mmean = 0
    M2mean = (8*(1 + exp(8*β)))/(3 + cosh(8*β))
    Mabsmean = (2*(2 + exp(8β)))/(3 + cosh(8*β))

    χ = (M2mean - Mabsmean^2)/(4*T) # magnetic susceptibility

    """ per particle """
    Emean /= 4
    Mmean /= 4
    Mabsmean /= 4

    expecvalues = [Emean, Cv, Mmean, χ, Mabsmean]

    if ssave
        savepath = string(datapath, "analytical_T=$(Float64(T))", ".jld")
        save(savepath, "expecvalues", expecvalues)
    end

    return expecvalues
end

function energy(lattice::Array, L::Real)
    """ Computes the energy in the current configuration
    """

    E = 0.0
    L = Int16(L)
    for j = 2:L+1, i = 2:L+1
        E -= lattice[i, j]*(lattice[i-1, j] +
                            lattice[i, j+1])
    end

    return E
end

function create_ghost_points(array)

    """ Function for adding ghost points to spin lattice """

    top = reshape(array[1,:], (1, size(array)[1]))
    bottom = reshape(array[end,:], (1, size(array)[1]))

    left = reshape(array[:,1], (size(array)[2], 1))
    right = reshape(array[:,end], (size(array)[2], 1))

    left = vcat(0, left, 0)
    right = vcat(0, right, 0)

    # contruct new lattice
    lattice = vcat(bottom, array, top)
    lattice = hcat(right, lattice, left)

    return lattice
end

function init(L::Real, T::Float64, ordered=false)
    """ Initialize energy, spin lattice, and magnetization
    Input:
    -----------
    L: number of spins in each dimension
    T: temperature
    ordered: true/false- whether to initialize with all spins = 1

    Output:
    -----------------
    - spin lattice with ghost points of size (L+2)x(L+2)
    - initial energy
    - initial magnetization
    """

    """ Spin lattice """
    spins = zeros(Int8, L, L)
    if ordered
        fill!(spins, 1)
    else
        rand!(spins, -1:2:1)
    end


    # create ghost points for periodic conditions
    lattice = create_ghost_points(spins)

    """ Energy """
    E = energy(lattice, L)

    """ Magnetization """
    M = sum(spins)

    return lattice, E, M
end

function get_filepath(T, N, L)
    filepath = string(datapath, "approx_T=$(Float64(T))_log10N=$(log10(N))_L=$(Float64(L)).jld")
    return filepath
end

function compute_quants(avgs, params, save=false)
    """
    Compute expectation values following N Monte Carlo cycles

    Input:
    -----
    avgs: array with [E, E^2, M, M^2, |M|]
    params: Temperature, Number of MC cycles used, lattice size as [T, N, L]
    save: save results in a JLD file

    Returns:
    --------
    - Mean energy expectation value ⟨E⟩
    - Specific heat capacity per particle Cv
    - Mean magnetization expecation value ⟨M⟩
    - Magnetic susceptibility per particle χ
    - Mean absolute value of magnetization ⟨|M|⟩
    """
    T::Float64, N, L::Float64 = params

    avgs /= N # average over N mc cycles

    Emean    = avgs[1]
    E2mean   = avgs[2]
    Mmean    = avgs[3]
    M2mean   = avgs[4]
    Mabsmean = avgs[5]

    """ compute heat capacity and susceptibility """
    Cv = (E2mean - Emean*Emean)/(L^2*T^2)
    Chi = (M2mean - Mabsmean*Mabsmean)/(L^2*T)

    # per spin
    Emean /= L^2
    Mmean /= L^2
    Mabsmean /= L^2

    """ save results if requested """
    if save

        filepath = get_filepath(T, N, L)

        save(filepath, "Emean", Emean, "Cv", Cv,
                       "Chi", Chi, "Mabs", Mabsmean)
    end

    return [Emean, Cv, Mmean, Chi, Mabsmean]
end


function montecarlo(T, N, L, ordered=false)

    """
    Perform N Monte Carlo simulations with Metropolis sampling

    Input:
    -------
    T:       temperature
    N:       number of cycles
    L:       number of spins in each dimension
    ordered: true/false, whether to start in an ordered initial configuration (all spins = 1)

    Returns:
    ---------
    - summed up E, E^2, M, M^2, |M|
    """

    avgs = zeros(Float64, 5)  # array for storing expectation valyues
    prob = zeros(Float64, 17) # array for storing probability ratios for each ΔE


    ΔEs = [-8., -4., 0., 4., 8.] # possible ΔE values
    prob[[1, 5, 9, 13, 17]] = [exp(-dE/T) for dE in ΔEs] # probability ratio for each ΔE

    prob = SVector{17}(prob) # static array for speed

    """ initialize spin matrix, energy and magnetization """
    lattice, E, M = init(L, T, ordered)

    Pe = [] # for counting energies, if required
    accepts = 0 # for counting no. of accepted configs
    LL = L*L

    xy = collect(2:L+1) # lattice indices

    """ Begin N Monte Carlo simulations """
    @inbounds for iters = 1:N
        """ randomly shuffle for each iteration such
            that we go over all spins but in random order """
        x = shuffle(xy)
        y = shuffle(xy)

        """ loop over over spins and flip spin """
        @inbounds for (ix, iy) in zip(x, y)
            ΔE = 2*lattice[iy, ix]*(lattice[iy, ix+1] +
                                  lattice[iy, ix-1] +
                                  lattice[iy-1, ix] +
                                  lattice[iy+1, ix])

            if (ΔE > 0)

                if rand() <= prob[ΔE + 9]
                    lattice[iy, ix] *= -1     # flip spin
                    M += 2*lattice[iy, ix]    # update magnetization
                    E += ΔE                 # update energy
                    accepts += 1

                    """ uncomment if we want to count energies"""
                    # if iters > Int64(1e4)
                    #     push!(Pe, E)
                    # end
                else
                    continue
                end # end inner if

            else """ accept if ΔE <= 0"""
                lattice[iy, ix] *= -1     # flip spin
                M += 2*lattice[iy, ix]    # update magnetization
                E += ΔE                 # update energy
                accepts += 1
            end # end outer if

        end # end metropolis


    """ Update expectation values """
    avgs[1] += E
    avgs[2] += E*E
    avgs[3] += M
    avgs[4] += M*M
    avgs[5] += abs(M)

    end # end Monte Carlo


    return avgs, accepts#, Pe
end
