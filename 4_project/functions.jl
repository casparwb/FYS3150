using Random, StaticArrays


function energy(lattice::Array, L::Real)
    """ Computes the energy in the current configuration """

    E = 0.0
    L = Int16(L)
    for j = 2:L+1, i = 2:L+1
        E -= lattice[i, j]*(lattice[i-1, j] +
                            lattice[i, j+1])

    end

    return E
end

function create_ghost_points(array)

    """ Function for adding ghost points to an array """

    top = reshape(array[1,:], (1, size(array)[1]))
    bottom = reshape(array[end,:], (1, size(array)[1]))

    left = reshape(array[:,1], (size(array)[2], 1))
    right = reshape(array[:,end], (size(array)[2], 1))

    left = vcat(0, left, 0)
    right = vcat(0, right, 0)

    lattice = vcat(bottom, array, top)
    lattice = hcat(right, lattice, left)

    return lattice
end

function init(L::Real, T::Float64, ordered=false)
    """ Initialize energy, spin lattice, and magnetization """

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

function compute_quants(avgs, params, absvar=false)
    """
    Compute expectation values

    Input:
    -----
    avgs: array with [E, E^2, M, M^2, |M|]
    params: Temperature, Number of MC cycles used, lattice size as [T, N, L]
    save: save results in a JLD file
    absvar: use |M| to compute variance in M

    Returns:
    --------
    - Mean energy expectation value ⟨E⟩
    - Specific heat capacity per particle Cv (var(E))
    - Mean magnetization expecation value ⟨M⟩
    - Magnetic susceptibility per particle χ (var(M))
    - Mean absolute value of magnetization ⟨|M|⟩
    """
    T::Float64, N, L::Float64 = params

    avgs /= N # average over N mc cycles

    Emean = avgs[1]
    E2mean = avgs[2]
    Mmean = avgs[3]
    M2mean = avgs[4]
    Mabsmean = avgs[5]

    Evar = (E2mean - Emean*Emean)/(L^2*T^2)
    if absvar
        Mvar = (M2mean - Mabsmean*Mabsmean)/(L^2*T)
    else
        Mvar = (M2mean - Mmean*Mmean)/(L^2*T)
    end


    # per spin
    Emean /= L^2
    Mmean /= L^2
    Mabsmean /= L^2

    # if save
    #
    #     filepath = get_filepath(T, N, L)
    #
    #     save(filepath, "Emean", Emean, "Evar", Evar,
    #                    "Mvar", Mvar, "Mabs", Mabsmean)
    # end

    return [Emean, Evar, Mmean, Mvar, Mabsmean]
end

function montecarlo(T, N, L, ordered=false)

    """
    Perform N Monte Carlo simulations with Metropolis sampling
    Input
    -------
    T:       temperature
    N:       number of cycles
    L:       number of spins in each dimension

    Returns:
    ---------
    - expectation value at every 100 time step
    """

    avgs = zeros(Float64, 5)
    Es = zeros(Float64, 17)

    ΔEs = [-8., -4., 0., 4., 8.]

    Es[[1, 5, 9, 13, 17]] = [exp(-E/T) for E in ΔEs]

    Es = SVector{17}(Es) # static array for speed

    """ initialize spin matrix, energy and magnetization """
    lattice, E, M = init(L, T, ordered)

    Pe = [] # for counting energies, if required
    accepts = 0 # for counting no. of accepted configs
    L2 = L*L
    """ Perform N Monte Carlo simulations """
    @inbounds for iters = 1:N

        """ Metropolis """
        @inbounds for n = 1:L2
            x = rand(2:L+1)
            y = rand(2:L+1)

            ΔE = 2*lattice[y, x]*(lattice[y, x+1] +
                                  lattice[y, x-1] +
                                  lattice[y-1, x] +
                                  lattice[y+1, x])



            if rand() <= Es[ΔE + 9]
                lattice[y, x] *= -1
                M += 2*lattice[y, x]
                E += ΔE
                accepts += 1

                """ uncomment if we want to count energies"""
                # if iters > Int64(1e4)
                #     push!(Pe, E)
                # end

            end # end if

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
