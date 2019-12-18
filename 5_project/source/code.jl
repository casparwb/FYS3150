""" Disease Modeling using the SIRS model """

using Revise
using PyPlot

include("functions.jl")

function f1(u, t, params)
    """ Simple model
    --------------------
    u: array with PDEs such that
        u[1] = S'(t)
        u[2] = I'(t)
        u[3] = R'(t)
    """
    a, b, c = params[1:3]

    N = sum(u)
    return [c*(u[3]) - a*u[1]*u[2]/N,
            a*u[1]*u[2]/N - b*u[2],
            b*u[2] - c*u[3]]
end

function f2(u, t, params)
    """ Model with vital dynamics
    -----------------------------
    u: array with PDEs such that
        u[1] = S'(t)
        u[2] = I'(t)
        u[3] = R'(t)
    """

    a, b, c, e, d, dI = params[1:6]

    N = sum(u)
    return [c*u[3] - a*u[1]*u[2]/N - d*u[1] + e*N,      # S'(t)
            a*u[1]*u[2]/N - b*u[2] - d*u[2] - dI*u[2],  # I'(t)
            b*u[2] - c*u[3] - d*u[3]]                   # R'(t)
end

function f3(u, t, params)
    """ Model with vital dynamics and seasonal variation
    --------------------------------------------------
    u: array with PDEs such that
        u[1] = S'(t)
        u[2] = I'(t)
        u[3] = R'(t)
    """

    a0, b, c, e, d, dI, A, ω = params[1:8]

    a = A*cos(ω*t) + a0
    N = sum(u)
    return [c*u[3] - a*u[1]*u[2]/N - d*u[1] + e*N,      # S'(t)
            a*u[1]*u[2]/N - b*u[2] - d*u[2] - dI*u[2],  # I'(t)
            b*u[2] - c*u[3] - d*u[3]]                   # R'(t)
end

function f4(u, t, params)
    """ Model with vaccination (all parameters)
    ----------------------------
    u: array with PDEs such that
        u[1] = S'(t)
        u[2] = I'(t)
        u[3] = R'(t)
    """
    a0, b, c, e, d, dI, A, ω, f = params

    a = A*cos(ω*t) + a0
    N = sum(u)
    return [c*u[3] - a*u[1]*u[2]/N - d*u[1] + e*N - f*u[1],      # S'(t)
            a*u[1]*u[2]/N - b*u[2] - d*u[2] - dI*u[2],           # I'(t)
            b*u[2] - c*u[3] - d*u[3] + f*u[1]]                   # R'(t)
end



function plotSIR_RK4(pops, filename, title)
    """ Plotting the populations computed with the RK4 method
    pops: arrays with each population, of shape (4, n, n) [population, u, t]"""

    pops = [pops[1] pops[2]; pops[3] pops[4]]

    t = pops[1, 1][2]

    labels = ["A" "B"; "C" "D"]

    fig, ax = subplots(2, 2, sharex=true)

    neqs = size(pops[1,1][1])[2] # number of equations


    for i = 1:2, j = 1:2
        pop = pops[i, j]
        S = pop[1][:,1]
        I = pop[1][:,2]
        R = pop[1][:,3]


        if i == 1 && j == 2
            ax[i,j].plot(t, S, label="S")
            ax[i,j].plot(t, I, label="I")
            ax[i,j].plot(t, R, label="R")
            ax[i,j].plot(t, (S + I + R), label="S + I + R")
            ax[i,j].set_title(labels[i, j])
        else
            ax[i,j].plot(t, S)
            ax[i,j].plot(t, I)
            ax[i,j].plot(t, R)
            ax[i,j].plot(t, (S + I + R))
            ax[i,j].set_title(labels[i, j])
        end


    end

    fig.text(0.5, 0.04, "Time", ha="center")
    fig.text(0.04, 0.5, "Counts", va="center", rotation="vertical")

    fig.legend(loc="center right")
    suptitle(title, y=1.006)
    figpath = "figures/" * filename * ".png"
    fig.savefig(figpath, dpi=500)
    gcf()

end

function plotSIR_MC(pops, t, filename, title)
    """ Plot the Monte Carlo results """
    A, B, C, D = pops

    fig, ax = subplots(2, 2, sharex=false)#, sharey=true)

    Na = Int(length(t)/4)
    S = A[1,:]
    I = A[2,:]
    R = A[3,:]

    ax[1,1].plot(t[1:Na], S[1:Na])
    ax[1,1].plot(t[1:Na], I[1:Na])
    ax[1,1].plot(t[1:Na], R[1:Na])
    ax[1,1].plot(t[1:Na], (S + I + R)[1:Na])
    ax[1,1].set_title("A")

    Nb = Int(length(t)/2)
    S = B[1,:]
    I = B[2,:]
    R = B[3,:]

    ax[1,2].plot(t[1:Nb], S[1:Nb], label="S")
    ax[1,2].plot(t[1:Nb], I[1:Nb], label="I")
    ax[1,2].plot(t[1:Nb], R[1:Nb], label="R")
    ax[1,2].plot(t[1:Nb], (S + I + R)[1:Nb], label="S + I + R")
    ax[1,2].set_title("B")

    Nc = Int(3*length(t)/4)
    S = C[1,:]
    I = C[2,:]
    R = C[3,:]

    ax[2,1].plot(t[1:Nc], S[1:Nc])
    ax[2,1].plot(t[1:Nc], I[1:Nc])
    ax[2,1].plot(t[1:Nc], R[1:Nc])
    ax[2,1].plot(t[1:Nc], (S + I + R)[1:Nc])
    ax[2,1].set_title("C")

    S = D[1,:]
    I = D[2,:]
    R = D[3,:]

    ax[2,2].plot(t, S)
    ax[2,2].plot(t, I)
    ax[2,2].plot(t, R)
    ax[2,2].plot(t, (S + I + R))
    ax[2,2].set_title("D")

    fig.text(0.5, 0.04, "Time", ha="center")
    fig.text(0.04, 0.5, "Counts", va="center", rotation="vertical")

    fig.legend(loc="center right")
    figpath = "figures/" * filename * ".png"
    suptitle(title, y=1.006)
    fig.savefig(figpath, dpi=500)
    gcf()

end

function simple_model(inits, params, T, n)
    """ Simulate the simple model without any added complexities
    params:[a, b, c, 0, 0, 0, 0, 0, 0]"""

    names  = ["A", "B", "C", "D"]

    bs = [params[i][2] for i = 1:4]
    a, c = params[1][1:2:3]

    exacts = [[b/a, (1 - b/a)/(1 + b/c), b/c*(1 - b/a)/(1 + b/c)] for b in bs]

    N = sum(inits)


    """ RK4 """
    #parameters for each population
    Aparams, Bparams, Cparams, Dparams = params

    A, t1 = @timed solve_RK4(f1, inits, T, n, Aparams)
    B, t2 = @timed solve_RK4(f1, inits, T, n, Bparams)
    C, t3 = @timed solve_RK4(f1, inits, T, n, Cparams)
    D, t4 = @timed solve_RK4(f1, inits, T, n, Dparams)

    t_rk4 = (t1 + t2 + t3 + t4)/4.

    title = "RK4 Simple Model"
    plotSIR_RK4([A, B, C, D], "RK4_simple_model", title)


    """ Monte Carlo """
    #inits = (S0, I0, N - (S0 + I0))

    M = 40000
    A = zeros(Float64, 3, M)
    B = zeros(Float64, 3, M)
    C = zeros(Float64, 3, M)
    D = zeros(Float64, 3, M)
    t = collect(1:M)
    pops = []

    for j = 1:M
        A[:,j] = MonteCarlo_SIR(j, T, inits, Aparams)[1][1:3]
        B[:,j] = MonteCarlo_SIR(j, T, inits, Bparams)[1][1:3]
        C[:,j] = MonteCarlo_SIR(j, T, inits, Cparams)[1][1:3]
        D[:,j] = MonteCarlo_SIR(j, T, inits, Dparams)[1][1:3]
    end

    title = "Monte Carlo Simple Model"
    plotSIR_MC([A, B, C, D], t, "MC_simple_model", title)

    expecsA, t1 = @timed MonteCarlo_SIR(M, T, inits, Aparams, true)[1]
    expecsB, t2 = @timed MonteCarlo_SIR(M, T, inits, Bparams, true)[1]
    expecsC, t3 = @timed MonteCarlo_SIR(M, T, inits, Cparams, true)[1]
    expecsD, t4 = @timed MonteCarlo_SIR(M, T, inits, Dparams, true)[1]

    t_mc = (t1 + t2 + t3 + t4)/4.

    expecs = [expecsA, expecsB, expecsC, expecsD]

    println("With $M Monte Carlo cycles (including 2000 burn in):")
    for (idx, nam) in enumerate(names)
        S, I, R, S2, I2, R2 = expecs[idx] ./ N
        varS = round(variance(S, S2), digits=2)
        varI = round(variance(I, I2), digits=2)
        varR = round(variance(R, R2), digits=2)

        println("$nam:")
        println("S = $(round(S, digits=2)) | varS = $varS")
        println("I = $(round(I, digits=2)) | varI = $varI")
        println("R = $(round(R, digits=2)) | varR = $varR")
        println("Exact [S, I, R] = $(exacts[idx])")
        println()
    end

    println("Times: \n RK4: $(t_rk4) | MC (M = $M): $(t_mc)")

    """
    julia> include("code.jl")
    With 40000 Monte Carlo cycles (including 2000 burn in):
    A:
    S = 0.25 | varS = 26.23
    I = 0.25 | varI = 25.1
    R = 0.5 | varR = 98.34
    Exact [S, I, R] = [0.25, 0.25, 0.5]

    B:
    S = 0.52 | varS = 112.56
    I = 0.09 | varI = 3.88
    R = 0.38 | varR = 60.36
    Exact [S, I, R] = [0.5, 0.1, 0.4]

    C:
    S = 0.96 | varS = 371.64
    I = 0.0 | varI = 0.02
    R = 0.04 | varR = 4.61
    Exact [S, I, R] = [0.75, 0.0357143, 0.214286]

    D:
    S = 0.96 | varS = 373.98
    I = 0.0 | varI = 0.01
    R = 0.03 | varR = 2.65
    Exact [S, I, R] = [1.0, 0.0, 0.0]

    Times:
     RK4: 0.006209899749999999 | MC (M = 40000): 0.003632725
     """
end

function vital_dynamics(inits, params, T, n)
    """ Study the effect of vital dynamics
    params: [a, b, c, e, d, dI, 0, 0, 0]"""

    """ RK4 """
    #parameters for each population
    Aparams, Bparams, Cparams, Dparams = params

    e, d, dI = Aparams[4:6]

    A = solve_RK4(f2, inits, T, n, Aparams)
    B = solve_RK4(f2, inits, T, n, Bparams)
    C = solve_RK4(f2, inits, T, n, Cparams)
    D = solve_RK4(f2, inits, T, n, Dparams)


    title = "RK4 w/ Vital Dynamics"
    plotSIR_RK4([A, B, C, D], "RK4_vital_dynamics_$dI", title)


    """ Monte Carlo """

    M = 40000
    A = zeros(Float64, 3, M)
    B = zeros(Float64, 3, M)
    C = zeros(Float64, 3, M)
    D = zeros(Float64, 3, M)
    t = collect(1:M)
    pops = []

    for j = 1:M
        A[:,j] = MonteCarlo_SIR(j, T, inits, Aparams)[1][1:3]
        B[:,j] = MonteCarlo_SIR(j, T, inits, Bparams)[1][1:3]
        C[:,j] = MonteCarlo_SIR(j, T, inits, Cparams)[1][1:3]
        D[:,j] = MonteCarlo_SIR(j, T, inits, Dparams)[1][1:3]
    end

    title = "Monte Carlo w/ Vital Dynamics \n e = $e, d = $d, dI = $dI"
    plotSIR_MC([A, B, C, D], t, "MC_vital_dynamics_$dI", title)


end

function seasonal_variation(inits, params, T, n)
    """ Study the effect of seasonal variations
    params: [a, b, c, d, dI, A0, ω, 0]"""


    #parameters for each population
    Aparams, Bparams, Cparams, Dparams = params

    e, d, dI, A0, ω = Aparams[4:8]

    A = solve_RK4(f3, inits, T, n, Aparams)
    B = solve_RK4(f3, inits, T, n, Bparams)
    C = solve_RK4(f3, inits, T, n, Cparams)
    D = solve_RK4(f3, inits, T, n, Dparams)


    title = "RK4 w/ seasonal variation \n e = $e, d = $d, dI = $dI, A0 = $A0, omega = $(round(ω, digits=2))\n"
    plotSIR_RK4([A, B, C, D], "RK4_seasonal_variation_A0$A0, omega$(round(ω, digits=2))", title)


    """ Monte Carlo """

    M = 40000
    A = zeros(Float64, 3, M)
    B = zeros(Float64, 3, M)
    C = zeros(Float64, 3, M)
    D = zeros(Float64, 3, M)
    t = collect(1:M)

    for j = 1:M
        A[:,j] = MonteCarlo_SIR(j, T, inits, Aparams)[1][1:3]
        B[:,j] = MonteCarlo_SIR(j, T, inits, Bparams)[1][1:3]
        C[:,j] = MonteCarlo_SIR(j, T, inits, Cparams)[1][1:3]
        D[:,j] = MonteCarlo_SIR(j, T, inits, Dparams)[1][1:3]
    end

    title = "Monte Carlo w/ seasonal variation \n e = $e, d = $d, dI = $dI, A0 = $A0, omega = $(round(ω, digits=2))\n"
    plotSIR_MC([A, B, C, D], t, "MC_seasonal_variation_A0$A0, omega$(round(ω, digits=2))", title)

end


function vaccination(inits, params, T, n)
    """ Study model with all added parameters """

    """ RK4 """
    #parameters for each population
    Aparams, Bparams, Cparams, Dparams = params

    e, d, dI, A0, ω, f = Aparams[4:end]

    A = solve_RK4(f4, inits, T, n, Aparams)
    B = solve_RK4(f4, inits, T, n, Bparams)
    C = solve_RK4(f4, inits, T, n, Cparams)
    D = solve_RK4(f4, inits, T, n, Dparams)

    title = "RK4 w/ vaccination \n e = $e, d = $d, dI = $dI, A0 = $A0, omega = $(round(ω, digits=2)), f=$f\n"
    plotSIR_RK4([A, B, C, D], "RK4_vaccination_dI$(dI)_A0$(A0)_w$(round(ω, digits=2))_f$f", title)


    """ Monte Carlo """

    M = 40000
    A = zeros(Float64, 3, M)
    B = zeros(Float64, 3, M)
    C = zeros(Float64, 3, M)
    D = zeros(Float64, 3, M)
    t = collect(1:M)

    for j = 1:M
        A[:,j] = MonteCarlo_SIR(j, T, inits, Aparams)[1][1:3]
        B[:,j] = MonteCarlo_SIR(j, T, inits, Bparams)[1][1:3]
        C[:,j] = MonteCarlo_SIR(j, T, inits, Cparams)[1][1:3]
        D[:,j] = MonteCarlo_SIR(j, T, inits, Dparams)[1][1:3]
    end

    title = "Monte Carlo w/ vaccination \n e = $e, d = $d, dI = $dI, A0 = $A0, omega = $(round(ω, digits=2)), f=$f\n"
    plotSIR_MC([A, B, C, D], t, "MC_vaccination_dI$(dI)_A0$(A0)_w$(round(ω, digits=2))_f$f", title)


end


function study_time(T=20)
    """ Study the CPU time and error difference  between the Runge Kutta
    method and the Monte Carlo implementation """
    N = 400
    S0 = 300
    I0 = 100
    R0 = 0

    inits = [S0, I0, R0]

    bs = [1, 2, 3, 4] # b-values
    a = 4
    c = 0.5

    """ use the simple model """
    Aparams = [a, bs[1], c, 0, 0, 0, 0, 0, 0]
    Bparams = [a, bs[2], c, 0, 0, 0, 0, 0, 0]
    Cparams = [a, bs[3], c, 0, 0, 0, 0, 0, 0]
    Dparams = [a, bs[4], c, 0, 0, 0, 0, 0, 0]

    params = [Aparams, Bparams, Cparams, Dparams]

    exacts = [[b/a, (1 - b/a)/(1 + b/c), (b/c)*(1 - b/a)/(1 + b/c)] for b in bs] # exact values of S, I, R in each population

    ns = [100, 500, 1000, 2500, 5000]
    Ms = [1000, 5000, 10000, 50000, 100000]

    """ arrays for storing CPU times and errors """
    times_RK = zeros(length(ns))
    times_MC = zeros(length(ns))

    errors_RK = zeros(3, 4, length(ns))
    errors_MC = zeros(3, 4, length(ns))

    for i = 1:length(ns)
        println("Starting i = $i/$(length(ns))")

        for j = 1:4 # each population
            RK4_result, tRK = @timed solve_RK4(f1, inits, T, ns[i], params[j])[1]
            MC_result, tMC = @timed MonteCarlo_SIR(Ms[i], T, inits, params[j], true)[1][1:3]

            RK4_result /= N # get fraction
            MC_result /= N

            times_RK[i] += tRK
            times_MC[i] += tMC

            for k = 1:3 # each group (S, I, R)
                errors_RK[k, j, i] = abs(((RK4_result[end, k] - exacts[j][k]) / exacts[j][k]))
                errors_MC[k, j, i] = abs(((MC_result[k] - exacts[j][k]) / exacts[j][k]))
            end # end k
        end # end j

        # average
        times_RK[i] /= 4
        times_MC[i] /= 4

    end # end i

    """ plot Monte Carlo results """
    fig = figure()
    subplot(221)
    plot(log10.(times_MC), log10.(errors_MC[1,1,:] ./ times_MC))
    plot(log10.(times_MC), log10.(errors_MC[2,1,:] ./ times_MC))
    plot(log10.(times_MC), log10.(errors_MC[3,1,:] ./ times_MC))
    title("A")

    subplot(222)

    plot(log10.(times_MC), log10.(errors_MC[1,2,:] ./ times_MC))
    plot(log10.(times_MC), log10.(errors_MC[2,2,:] ./ times_MC))
    plot(log10.(times_MC), log10.(errors_MC[3,2,:] ./ times_MC))
    title("B")

    subplot(223)

    plot(log10.(times_MC), log10.(errors_MC[1,3,:] ./ times_MC))
    plot(log10.(times_MC), log10.(errors_MC[2,3,:] ./ times_MC))
    plot(log10.(times_MC), log10.(errors_MC[3,3,:] ./ times_MC))
    xlabel("C")

    subplot(224)

    plot(log10.(times_MC), log10.(errors_MC[1,4,:] ./ times_MC), label="S")
    plot(log10.(times_MC), log10.(errors_MC[2,4,:] ./ times_MC), label="I")
    plot(log10.(times_MC), log10.(errors_MC[3,4,:] ./ times_MC), label="R")
    xlabel("D")

    fig.text(0.5, 0.04, "log10 CPU time", ha="center")
    fig.text(0.04, 0.5, "log10 Error / CPU time [s^-1]", va="center", rotation="vertical")

    suptitle("Monte Carlo error per CPU time as function of CPU time")
    legend()
    savefig("figures/MC_error_.png", dpi=500)

    """ plot Runge Kutta results """
    fig = figure()
    subplot(221)
    plot(log10.(times_RK), log10.(errors_RK[1,1,:] ./ times_RK))
    plot(log10.(times_RK), log10.(errors_RK[2,1,:] ./ times_RK))
    plot(log10.(times_RK), log10.(errors_RK[3,1,:] ./ times_RK))
    title("A")


    subplot(222)
    plot(log10.(times_RK), log10.(errors_RK[1,2,:] ./ times_RK))
    plot(log10.(times_RK), log10.(errors_RK[2,2,:] ./ times_RK))
    plot(log10.(times_RK), log10.(errors_RK[3,2,:] ./ times_RK))

    title("B")

    subplot(223)
    plot(log10.(times_RK), log10.(errors_RK[1,3,:] ./ times_RK))
    plot(log10.(times_RK), log10.(errors_RK[2,3,:] ./ times_RK))
    plot(log10.(times_RK), log10.(errors_RK[3,3,:] ./ times_RK))
    xlabel("C")

    subplot(224)
    plot(log10.(times_RK), log10.(errors_RK[1,4,:] ./ times_RK), label="S")
    plot(log10.(times_RK), log10.(errors_RK[2,4,:] ./ times_RK), label="I")
    plot(log10.(times_RK), log10.(errors_RK[3,4,:] ./ times_RK), label="R")
    xlabel("D")

    fig.text(0.5, 0.04, "log10 CPU time", ha="center")
    fig.text(0.04, 0.5, "log10 Error / CPU time [s^-1]", va="center", rotation="vertical")

    suptitle("Runge Kutta error per CPU time as function of CPU time")
    legend()
    savefig("figures/RK_error_.png", dpi=500)


    gcf()
    close(fig)

end

function main(T, n)
    """ Main program used to produce the various results in the report """

    """ Constants """
    N = 400
    S0 = 300
    I0 = 100
    R0 = 0
    inits = [S0, I0, R0]

    bs = [1, 2, 3, 4] # b-values
    a = 4
    c = 0.5

    e = 0.02        # birth rate
    d = 0.01        # natural death rate
    ω = 2*π/(T/2)   # seasonal variation frequency

    """ parameters to vary """
    dIs = [0, 0.5, 1, 2]
    A0s = [1, 2, 3, 4]
    fs = [0.2, 0.4, 0.8, 1]

    """ simple model """

    Aparams = [a, bs[1], c, 0, 0, 0, 0, 0, 0]
    Bparams = [a, bs[2], c, 0, 0, 0, 0, 0, 0]
    Cparams = [a, bs[3], c, 0, 0, 0, 0, 0, 0]
    Dparams = [a, bs[4], c, 0, 0, 0, 0, 0, 0]
    simpleparams = [Aparams, Bparams, Cparams, Dparams]

    simple_model(inits, simpleparams, T, n)


    """ with vital dynamics """

    Aparams = [a, bs[1], c, e, d, dIs[2], 0, 0, 0]
    Bparams = [a, bs[2], c, e, d, dIs[2], 0, 0, 0]
    Cparams = [a, bs[3], c, e, d, dIs[2], 0, 0, 0]
    Dparams = [a, bs[4], c, e, d, dIs[2], 0, 0, 0]

    vital_params = [Aparams, Bparams, Cparams, Dparams]

    vital_dynamics(inits, vital_params, T, n)


    """ with seasonal variation """

    Aparams = [a, bs[1], c, e, d, dIs[2], A0s[1], ω, 0]
    Bparams = [a, bs[2], c, e, d, dIs[2], A0s[1], ω, 0]
    Cparams = [a, bs[3], c, e, d, dIs[2], A0s[1], ω, 0]
    Dparams = [a, bs[4], c, e, d, dIs[2], A0s[1], ω, 0]

    seasonal_params = [Aparams, Bparams, Cparams, Dparams]

    seasonal_variation(inits, seasonal_params, T, n)


    """ with vaccination """

    Aparams = [a, bs[1], c, e, d, dIs[2], A0s[1], ω, fs[3]]
    Bparams = [a, bs[2], c, e, d, dIs[2], A0s[1], ω, fs[3]]
    Cparams = [a, bs[3], c, e, d, dIs[2], A0s[1], ω, fs[3]]
    Dparams = [a, bs[4], c, e, d, dIs[2], A0s[1], ω, fs[3]]

    all_params = [Aparams, Bparams, Cparams, Dparams]

    vaccination(inits, all_params, T, n)


end
