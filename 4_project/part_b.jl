using Plots, LinearAlgebra, JLD, Random, StaticArrays, SharedArrays, Test, Revise
using PyCall
#using PlotThemes; theme(:juno, titlecolor="white")
const figpath = "figures/"
const datapath = "data/"

include("functions.jl")

function most_likely_state()
    """ When is most likely state reached? """
    L = 20
    T1 = 1.0
    T2 = 2.4

    cycles = collect(1:20:20000)
    expecs_ord = zeros(length(cycles), 2)
    expecs_unord = copy(expecs_ord)
    accepts = zeros(2, length(cycles), 2)
    for (Ti, Temp) in enumerate([T1, T2])
        println("Starting T=$Temp")
        for i = 1:length(cycles)
            avgs_ord, accepts[Ti, i, 1] = montecarlo(Temp, cycles[i], L, true)
            expecs_ord[i,:] = compute_quants(avgs_ord, [Temp, cycles[i], L])[1:3:4]

            avgs_unord, accepts[Ti, i, 2] = montecarlo(Temp, cycles[i], L, false)
            expecs_unord[i,:] = compute_quants(avgs_unord, [Temp, cycles[i], L])[1:3:4]

            if i % 10 == 0
                println("$i/$(length(cycles)) cycles = $(cycles[i]) finished")
            end
        end
        # plot energy

        # plot(cycles, expecs_ord[:,1], label="Ordered", legend=:bottomright)
        # plot!(cycles, expecs_unord[:,1], label="Unordered")
        # xlabel!("MC cycles")
        # ylabel!("Energy")
        # title!("T = $Temp")

        p1 = plot(cycles, expecs_ord[:,1], label="Ordered")
        p2 = plot(cycles, expecs_unord[:,1], label="Unordered", xlabel="MC cycles")
        plot(p1, p2,
             layout=(2, 1),
             ylabel="E",
             dpi=500)
        savefig("figures/new_energy_T=$(Float64(Temp)).png")

        # plot magnetization

        # plot(cycles, expecs_ord[:,end], label="Ordered")
        # plot!(cycles, expecs_unord[:,end], label="Unordered")
        # xlabel!("MC cycles")
        # ylabel!("Absolute Magnetization")
        # title!("T = $Temp")


        p3 = plot(cycles, expecs_ord[:,2], label="Ordered")
        p4 = plot(cycles, expecs_unord[:,2], label="Unordered", xlabel="MC cycles")
        plot(p3, p4,
             layout=(2, 1),
             ylabel="M",
             dpi=500)
        savefig("figures/new_absmag_T=$(Float64(Temp))")
    end

    accepts[1,:,1] ./= cycles
    accepts[1,:,2] ./= cycles
    accepts[2,:,1] ./= cycles
    accepts[2,:,2] ./= cycles

    p1 = plot(cycles, accepts[1,:,1], label="Ordered")
    plot!(cycles, accepts[1,:,2], label="Unordered")
    title!("T = $T1")

    p2 = plot(cycles, accepts[2,:,1], label="Ordered")
    plot!(cycles, accepts[2,:,2], label="Unordered")
    title!("T = $T2")

    plot(p1, p2, xlabel="MC cycles", ylabel="No. of accepted configs.")
    savefig("figures/new_counts.png")



end


function prob_dist(N)

    """
    Plot the probability distribution
    for the energy for T=1.0 and T=2.4
    """
    T1 = 1.0
    T2 = 2.4
    L = 20

    for Temp in [T1, T2]
        avgs, Pe = montecarlo(Temp, N, L)
        Pe /= L^2

        expecvals = compute_quants(avgs, [Temp, N, L])
        histogram(Pe,
                  normalize=:pdf,
                  legend=false,
                  dpi=500)
        if Temp == T1
            plot!([expecvals[1], expecvals[1]], [0, 600], linewidth=5)
        elseif Temp == T2
            plot!([expecvals[1], expecvals[1]], [0, 20], linewidth=5)
        end

        ylabel!("Probability density")
        xlabel!("Energy")
        title!("P(E, T=$Temp)")
        savefig("figures/new_prob_dist_T=$Temp.png")
    end
end


function main()
    """ Main program """

    # """ Comparing analytical values """
    # N = [10^i for i = 1:7]
    # analyticals = load("data/analytical_T=1.0.jld", "expecvalues") # [Emean, Evar, Mmean, Mvar, Mabs]
    # for i in 1:length(N)
    #     avgs = montecarlo(1.0, N[i], 2, true)
    #     expecs = compute_quants(avgs, [1.0, N[i], 2])
    #
    #     println()
    #     println("$expecs")
    # end
    # println("$analyticals")

    #prob_dist(N=1000000)
    most_likely_state()

    """ timing serial and comparing with parallel """
    # Ls = [40, 60]
    # temps = collect(2.0:0.05:2.3)
    #
    # times = zeros(len(temps), 2)
    # for (i, T) in enumerate(temps)
    #     for (j, L) in enumerate(Ls)
    #         avgs. times[i,1] = @timed montecarlo(T, N, L)
    #     end
    # end

    # times_serial = load("data/serial_times_new.jld", "times")
    # times_parallel = load("data/part_e.jld", "times")
    #
    # times_4060_serial = sum(times_serial, dims=1)
    # times_4060_parallel = [maximum(times_parallel[1:4,i]) for i = 1:2]
    # times_4060_parallel += [maximum(times_parallel[5:end,i]) for i = 1:2]
    #
    # for i = 1:2
    #     time_s = times_4060_serial[i]
    #     time_p = times_4060_parallel[i]
    #     speedup = time_s/time_p
    #
    #     println("Time serial: $(round(time_s)) | Time parallel: $(round(time_p))
    #             | Speedup: $(round(time_s/time_p))")
    # end

    """ studying phase transitions """
    # expecvals = load("data/part_e.jld", "expecvals")
    # """
    # - shape (7, 4, 4)
    # - dim 1: Temperatures 2.0:0.05:2.3
    # - dim 2: Spins L [40, 60, 80, 100]
    # - dim 3: E, Cv, χ, |M|
    # """
    # Temps = collect(2.0:0.05:2.3) # temperatues tested for
    # Ls = [40, 60, 80, 100]
    # Cvs = [expecvals[:,i,2] for i = 1:4]
    # Chis = [expecvals[:,i,3] for i = 1:4]
    #
    # np = pyimport("numpy")
    #
    # Temps_new = collect(2.0:0.001:2.3)
    #
    # """ fit Cv and χ to a polynomial for each L """
    # Cvs_new = []
    # Chis_new = []
    #
    # for i = 1:4
    #     polyCv = np.polyfit(Temps, Cvs[i], deg=4)
    #     polyChi = np.polyfit(Temps, Chis[i], deg=4)
    #
    #     new_Cv = np.polyval(polyCv, Temps_new)
    #     new_Chi = np.polyval(polyChi, Temps_new)
    #     push!(Cvs_new, new_Cv)
    #     push!(Chis_new, new_Chi)
    # end
    # # scatter(Temps, Cvs[1], label="Computed Values", dpi=500, legend=:bottomright)
    # # plot!(Temps_new, Cvs_new[1], label="Fitted values")
    # # xlabel!("T [kT/J]"); ylabel!("Specific Heat Capacity")
    # # title!("Fitted Cv values with 4th deg. polynomial")
    # # savefig("figures/fitted.png")
    #
    #
    # Tcrits_Cv = [Temps_new[argmax(Cvs_new[i])] for i = 1:4]
    # Tcrits_Chi = [Temps_new[argmax(Chis_new[i])] for i = 1:4]
    # Tcrits = [np.mean([Tcv, Tchi]) for (Tcv, Tchi) in zip(Tcrits_Cv, Tcrits_Chi)]
    #
    # Ls_inv = 1.0 ./ Ls
    # Ls_inv_new = collect(Ls_inv[end]:0.001:Ls_inv[1])
    #
    # poly = np.polyfit(Ls_inv, Tcrits_Cv, deg=1)
    # intercept = poly[end]
    # println(intercept)
    #
    # Tcrits_new = np.polyval(poly, Ls_inv_new)
    # scatter(Ls_inv, Tcrits_Cv, label="Computed Tc", dpi=500)
    # plot!(Ls_inv_new, Tcrits_new, label="Fitted Tc")
    # xlabel!("L^-1"); ylabel!("Critical Temperature Tc")
    # savefig("figures/crit_temp.png")

end
