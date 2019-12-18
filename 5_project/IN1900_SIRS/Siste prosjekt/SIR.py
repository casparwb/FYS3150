# Simulating the spreading of a disease by a SIR model

import numpy as np

def f_user(u, t, beta):
    """Defining the system of ODE's"""
    return [-beta*u[0]*u[1],       # S(t)
            beta*u[0]*u[1]-v*u[1], # I(t)
            u[1]*v]                # R(t)

def solve(f_user, U0, T, n):
    """Solving the system of differential equations
    using the RungeKutta4-method"""
    def f(u, t):
        return np.asarray(f_user(u, t, beta)) # Returning the pre-defined system of ODE's as an array
    t = np.zeros(n+1)                         # t-values
    U0 = np.asarray(U0)                       # Initial conditions
    neq = U0.size                             # Amount of columns
    u = np.zeros((n+1, neq))                  # An array with n+1-length and neq amount of coulmns (number of ODE's)
    t[0] = 0                                  # Initial t-value
    u[0] = U0                                 # Initial conditions
    tol = 1e-9
    for k in range(n):
        K1 = dt*f(u[k], t[k])
        K2 = dt*f(u[k]+0.5*K1, t[k]+0.5*dt)
        K3 = dt*f(u[k]+0.5*K2, t[k]+0.5*dt)
        K4 = dt*f(u[k]+K3, t[k]+dt)
        u[k+1] = u[k] + (1/6)*(K1+2*K2+2*K3+K4)
        t[k+1] = dt+t[k]
        success = ((u[k][0]+u[k][1]+u[k][2]) - (np.sum(U0))) < tol
        assert success, "Bug in SIR"
    return u, t


def plot():
    """Visualizing the differential equations"""
    import matplotlib.pyplot as plt
    u, t = solve(f_user, U0, T, n)
    S, I, R = u[:,0], u[:,1], u[:,2]
    plt.plot(t, S, label = "Susceptibles")
    plt.plot(t, I, label = "Infected")
    plt.plot(t, R, label = "Recovered")
    plt.legend(fontsize = 30)
    plt.xlabel("Days", fontsize = 20); plt.ylabel("People", fontsize = 20)
    plt.show()



N = 1500          # Total population
S0 = 1500         # Initial value of susceptibles
I0 = 1            # Initial value of infected
R0 = 0            # Initial value of recovered
U0 = [S0, I0, R0] # Initial values
v = 0.1
dt = 0.5          # Time step
T = 60
n = int(T/dt)     # Number of steps
beta = 0.0005
plot()
beta = 0.0001
plot()

"""As opposed to when beta equals 0.0005, when beta is 0.0001, the amount of infected
never goes above a certain low value, and all the ones who get infected, end up recovering.
This is hard to visualize for only 60 days, but when we do a year, we see how this
is true."""

"""
run SIR.py
-Output is a plot-
"""
