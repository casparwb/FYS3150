function solve_RK4(f, u0, T, n, params)
    """ Solves the coupled PDE using RungeKutta 4

    u0: initial conditions (population)
    T: maximum time
    n: number of time steps
    params: array with relevant parameters
    """

    t = zeros(Float64, n+1)

    """ check if single or coupled equations """
    if isa(u0, Tuple) || isa(u0, Array)
        n_eq = size(u0)[1]              # number of equations
        u = zeros(Float64, n+1, n_eq)   # array for storing solutions at each time step
        u[1,:] = u0
    else
        n_eq = 1
        u = zeros(Float64, n+1, n_eq)
        u[1] = u0
    end

    t[1] = 0.

    """ solve using the RK4 method """
    dt = T/n
    for k = 1:n
        K1 = dt*f(u[k,:], t[k], params)
        K2 = dt*f(u[k,:]+0.5*K1, t[k]+0.5*dt, params)
        K3 = dt*f(u[k,:]+0.5*K2, t[k]+0.5*dt, params)
        K4 = dt*f(u[k,:]+K3, t[k]+dt, params)
        u[k+1,:] = u[k,:] + (1/6)*(K1+2*K2+2*K3+K4)
        t[k+1] = dt+t[k]
    end

    """ return solutions and time """
    return u, t

end


function MonteCarlo_SIR(M, T, inits, params, burn=false)
    """ Simulate the model with with M Monte Carlo cycles

    input:
    ----------------
    M: number of MC cycles
    inits: initial condition (S0, I0, R0)
    params: relevant parameters (a, b, c, e, d, dI, A0, ω, f)
    burn: true/false, burn in the first 2000 iterations to produce more accurate results

    output:
    ----------------
    - average values of S, I, R, S^2, I^2, R^2 over the M cycles
    - total simulation time (sum of Δt)
    """

    # initialize parameters and initial conditions
    a0, b, c, e, d, dI, A0, ω, f = params
    S0, I0, R0 = inits

    expecs = zeros(2*length(inits))     # for storing expectation values

    N = sum(inits)                      # initial total population count

    # initialize time, population counts and rate of transmission
    t = 0.
    S, I, R = inits
    a = A0 + a0

    dt = T/M # time step for the variational rate of transmission

    """ burn in the first 2000 iterations if true """
    if burn

        for m = 1:2000

            Δt = minimum([4.0/(a*N), 1.0/(b*N), 1.0/(c*N), 1.0/(e*N), 1.0/(d*N), 1.0/(dI*I)]) # time step


            """ probabilities for each move """
            P_S_I = (a*S*I)/N*Δt
            P_I_R = b*I*Δt
            P_R_S = c*R*Δt

            "probabilites for birth and death"
            P_birth = e*N*Δt
            P_death_S = d*S*Δt
            P_death_I = d*I*Δt
            P_death_R = d*R*Δt
            P_deathI = dI*I*Δt

            """ check P(S → I) """
            if P_S_I >= rand()
                S -= 1
                I += 1
            end

            """ check P(I → R) """
            if P_I_R >= rand()
                I -= 1
                R += 1
            end

            """ check P(R → S) """
            if P_R_S >= rand()
                R -= 1
                S += 1
            end

            """ check birth and death probabilites """
            P_birth   >= rand() ? S += 1 : nothing
            P_death_S >= rand() ? S -= 1 : nothing
            P_death_I >= rand() ? I -= 1 : nothing
            P_death_R >= rand() ? R -= 1 : nothing
            P_deathI  >= rand() ? R -= 1 : nothing

            """ update total population count """
            N = S + I + R

            a = A0*cos(ω*(Δt*m)) + a0

        end
    end

    """ continue and save expectation values """
    @inbounds for m = 2:M+1

        Δt = minimum([4.0/(a*N), 1.0/(b*N), 1.0/(c*N), 1.0/(e*N), 1.0/(d*N), 1.0/(dI*N), 1.0/(f*N)]) # time step


        """ probabilities for each move """
        P_S_I = (a*S*I/N)*Δt
        P_I_R = b*I*Δt
        P_R_S = c*R*Δt

        """ probability for move S → R from vaccination """
        P_S_R = f*S*Δt

        """ probabilites for birth and death """
        P_birth   = e* N*Δt # birth
        P_death_S = d* S*Δt # natural death in population S
        P_death_I = d* I*Δt # natural death in population I
        P_death_R = d* R*Δt # natural death in population R
        P_deathI  = dI*I*Δt # death in I as a result of the disease


        """ check P(S → I) """
        if P_S_I > rand()
            S -= 1
            I += 1
        end

        """ check P(S → R) """
        if P_S_R > rand()
            S -= 1
            R += 1
        end

        """ check P(I → R) """
        if P_I_R > rand()
            I -= 1
            R += 1
        end

        """ check P(R → S) """
        if P_R_S > rand()
            R -= 1
            S += 1
        end

        """ check birth and death probabilites """
        P_birth   > rand() ? S += 1 : nothing
        P_death_S > rand() ? S -= 1 : nothing
        P_death_I > rand() ? I -= 1 : nothing
        P_death_R > rand() ? R -= 1 : nothing
        P_deathI  > rand() ? I -= 1 : nothing


        """ update total population count """
        N = S + I + R

        """ update time and a (if variational) """
        t += Δt
        a = A0*cos(ω*(M*dt)) + a0

        """ skip contributions if something has gone wrong (happened some times)""" 
        if S < 0
            expecs[1] += 0
            expecs[2] += I
            expecs[3] += R
            expecs[4] += 0
            expecs[5] += I*I
            expecs[6] += R*R
        elseif I < 0
            expecs[1] += S
            expecs[2] += 0
            expecs[3] += R
            expecs[4] += S*S
            expecs[5] += 0
            expecs[6] += R*R
        elseif R < 0
            expecs[1] += S
            expecs[2] += I
            expecs[3] += 0
            expecs[4] += S*S
            expecs[5] += I*I
            expecs[6] += 0
        else
            expecs[1] += S
            expecs[2] += I
            expecs[3] += R
            expecs[4] += S*S
            expecs[5] += I*I
            expecs[6] += R*R
        end

    end # end m

    expecs /= M # average over all cycles

    return expecs, t
end

function variance(X, X2)
    """ Return variance of the quantity X given ⟨X⟩ and ⟨X^2⟩ """
    return X2 - X*X
end
