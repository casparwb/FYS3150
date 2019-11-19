using Plots, JLD
@everywhere using SharedArrays
@everywhere include("functions.jl")


function parallel(N)
    """
    Parallel implementation of phase transition study

    N: number of Monte Carlo cycles

    - saves resulting expectation values for each temperature, for each lattice size
      in a JLD files
    - also saves the times taken
    """

    # broadcast necessary data to all workers
    @everywhere begin
        Temps = collect(2.0:0.05:2.3) # temperature range to study
        Ls = [40, 60, 80, 100]        # lattice sizes

        expecvals = SharedArray{Float64}((length(Temps), length(Ls), 4))
        times = SharedArray{Float64}((length(Temps), length(Ls)))
        ntemps = length(Temps)
    end

    # begin parallelized (distributed) temperature loop
    @sync @distributed for i = 1:ntemps
        for j = 1:4
            avgs, times[i, j] = @timed montecarlo(Temps[i], N, Ls[j])
            expecvals[i, j, :] = compute_quants(avgs, [Temps[i], N, Ls[j]])
        end
    end

    # extract
    E = expecvals[:,:,1]
    Cv = expecvals[:,:,2]
    χ = expecvals[:,:,3]
    Mabs = expecvals[:,:,5]

    """ plot result and save figure """
    p1 = plot(legend=false); p2 = plot(legend=:topright); p3 = plot(legend=false); p4 = plot(legend=false)
    p = [p1, p2, p3, p4]
    labels = ["E", "Cv", "Chi", "Mabs"]
    gr()
    for i = 1:4
        p[i] = plot()
        for j = 1:4
            if i == 2
                plot!(Temps, expecvals[:,j,i], ylabel=labels[i], label="L = $(Ls[j])")
            else
                plot!(Temps, expecvals[:,j,i], ylabel=labels[i], legend=false)
            end
        end
    end

    plot(p[1], p[2], p[3], p[4], layout=(2, 2), dpi=500, legendfontsize=5)
    xlabel!("T [kT/J]")
    savefig("figures/part_f.png")


    # extract results and save
    expecvalss = sdata(expecvals)
    timess = sdata(times)

    save("data/part_e.jld", "expecvals", expecvalss, "times", timess)

end

# n = 100 2.664956 seconds (3.78 M allocations: 80.478 MiB, 0.92% gc time)
# n = 1000 1.988140 seconds (3.62 M allocations: 70.801 MiB, 0.97% gc time)
# n = 10000 15.712823 seconds (3.62 M allocations: 70.714 MiB, 0.09% gc time)
# n = 100000 149.865587 seconds (3.56 M allocations: 69.658 MiB, 0.01% gc time)
