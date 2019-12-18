#Simulating human-zombie interaction

"""Defining the groups"""

# S: susceptible humans who can become zombies
# I: Infected humans, being bitten by zombies
# Z: zombies
# R: removed individuels, conquered by zombies or dead humans

"""Defining the consequential ODE's to solve"""

# S' = sigma - beta*S*Z-delta_S*S
# I' = beta*S*Z-rho*I-delta_I*I
# Z' = rho*I-alpha*S*Z
# R' = delta_S*S+delta_I*I+alpha*S*Z


"""
Defining the variables:

Sigma - number of new humans brought into the zombified area per unit Time

beta - the propability that a theoritcally possible human-zombie pair actually
meets physically, during a unit time interval, with the result that the human
is infected

delta S(dS) - the propability that a susceptible human is killed or dies,
in a unit time interval

delta I(dI) - the propability that an infected human is killed or dies,
in a unit time interval

rho - the propability that an infected human is turned into a zombie,
during a unit time interval

alpha - the propability that, during a unit time interval,
a theoritcally possible human-zombie pair fights and the human kills
the zombie
"""

import numpy as np
import ODESolver

class ProblemSIZR(object):
    def __init__(self, sigma, beta, dS, dI, rho, alpha, S0, I0, Z0, R0, T):
        self.sigma, self.beta, self.dS, self.dI, self.rho, \
        self.alpha, self.S0, self.I0, self.Z0, self.R0, self.T = \
        sigma, beta, dS, dI, rho, alpha, S0, I0, Z0, R0, T

        """Allowing the variables to be functions of time"""

        if isinstance(sigma, (float,int)):
            self.sigma = lambda t: sigma
        elif callable(sigma):
            self.sigma = sigma

        if isinstance(beta, (float, int)):
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta

        if isinstance(dS, (float,int)):
            self.dS = lambda t: dS
        elif callable(dS):
            self.dS = dS

        if isinstance(dI, (float, int)):
            self.dI = lambda t: dI
        elif callable(dI):
            self.dI = dI

        if isinstance(rho, (float,int)):
            self.rho = lambda t: rho
        elif callable(rho):
            self.rho = rho

        if isinstance(alpha, (float, int)):
            self.alpha = lambda t: alpha
        elif callable(alpha):
            self.alpha = alpha

    def __call__(self, u, t):
        """The system of ODE's to solve"""
        S, I, Z, R = u
        return [self.sigma(t) - self.beta(t)*S*Z - self.dS(t)*S,  # S'(t)
                self.beta(t)*S*Z - self.rho(t)*I - self.dI(t)*I,  # I'(t)
                self.rho(t)*I - self.alpha(t)*S*Z,                # Z'(t)
                self.dS(t)*S + self.dI(t)*I + self.alpha(t)*S*Z]  # R'(t)

class SolverSIZR(object):
    """Solving the system of ODE's using the 4th order Runge-Kutta method
    from the ODESolver"""
    def __init__(self, problem, dt):
        self.problem, self.dt = problem, dt

    def solve(self, method=ODESolver.RungeKutta4):
        self.solver = method(self.problem)
        ic = [self.problem.S0, self.problem.I0, self.problem.Z0, self.problem.R0]
        self.solver.set_initial_condition(ic)
        n = int(round(self.problem.T/self.dt))
        t = np.linspace(0, self.problem.T, n+1)
        u, self.t = self.solver.solve(t)
        self.S, self.I, self.Z, self.R = u[:,0], u[:,1], u[:,2], u[:,3]
        return u, t

    def plot(self, title):
        """Visualizing the SIZR-system"""
        import matplotlib.pyplot as plt
        u, t = self.solve()
        S = u[:,0]; I = u[:,1]; Z = u[:,2]; R = u[:,3]
        plt.plot(t, S, label = "Susceptibles")
        plt.plot(t, I, label = "Infected")
        plt.plot(t, Z, label = "Zombies")
        plt.plot(t, R, label = "Removed")
        plt.legend(fontsize = 30)
        plt.xlabel("Hours", fontsize = 20); plt.ylabel("People", fontsize = 20)
        plt.title(title, fontsize = 30)
        plt.show()

beta = 0.0012; alpha = 0.0016; dI = 0.014; sigma = 2; rho = 1; dS = 0;
S0 = 10; Z0 = 100; I0 = 0; R0 = 0
T = 24

problem = ProblemSIZR(sigma, beta, dS, dI, rho, alpha, S0, Z0, I0, R0, T)
solver = SolverSIZR(problem, dt=0.25)
solver.plot("Human-zombie interaction: Hysterial phase")
