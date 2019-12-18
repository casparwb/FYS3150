"""Implementing the SIR-model in a class"""

import ODESolver
import numpy as np

class ProblemSIR(object):
    """Defining the problem of ODE's to be solves"""
    def __init__(self, nu, beta, S0, I0, R0, T):
        self.S0, self.I0, self.R0, self.T = S0, I0, R0, T

        """Allowing variables to be functions:"""

        if isinstance(nu, (float,int)):
            self.nu = lambda t: nu
        elif callable(nu):
            self.nu = nu

        if isinstance(beta, (float, int)):
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta

    def __call__(self, u, t):
        """Right-hand side function of the ODE system"""
        S, I, R = u
        return [-self.beta(t)*S*I,             # S'(t)
                self.beta(t)*S*I-self.nu(t)*I, # I'(t)
                self.nu(t)*I]                  # R'(t)

class SolverSIR(object):
    """Solving the system of ODE's"""
    def __init__(self, problem, dt):
        self.problem, self.dt = problem, dt

    def solve(self, method=ODESolver.RungeKutta4):
        """Using the pre-written ODESolver to solve the ODEs"""
        self.solver = method(self.problem)                       # An instance of the solver is created
        ic = [self.problem.S0, self.problem.I0, self.problem.R0] # The initial conditions are defined as a list
        self.solver.set_initial_condition(ic)                    # The list of initial conditions are set in the ODESolver
        n = int(round(self.problem.T/self.dt))                   # n is calculated
        t = np.linspace(0, self.problem.T, n+1)                  # t is calculated
        u, self.t = self.solver.solve(t)                         # u and t is computed and extracted
        self.S, self.I, self.R = u[:,0], u[:,1], u[:,2]
        return u, t

    def plot(self, title):
        import matplotlib.pyplot as plt
        u, t = self.solve()                    # We compute and extract u and t using the solver in the same class
        S = u[:,0]; I = u[:,1]; R = u[:,2]     # Extracting the different ODE's to plot
        plt.plot(t, S, label = "Susceptibles")
        plt.plot(t, I, label = "Infected")
        plt.plot(t, R, label = "Recovered")
        plt.legend(fontsize = 30)
        plt.title(title, fontsize=30)
        plt.xlabel("Days", fontsize = 20); plt.ylabel("People", fontsize = 20)
        plt.show()

problem = ProblemSIR(beta=lambda t: 0.0005 if t <= 12 else 0.0001, nu = 0.1, S0 = 1500, I0=1, R0=0, T=60)
solver = SolverSIR(problem, dt=0.5)
u = solver.solve()[0]
max_inf = max(u[:,1])

#Instancing with a constant beta = 0.0005:

problem2 = problem = ProblemSIR(beta=0.0005, nu = 0.1, S0 = 1500, I0=1, R0=0, T=60)
solver2 = SolverSIR(problem2, dt=0.5)
u2 = solver2.solve()[0]
max_inf2 = max(u2[:,1])

print("""With beta varying between 0.0005 and 0.0001, the maximum number of infected is %g.

With beta at a constant 0.0005, the maximum number of infected is %g""" %(max_inf, max_inf2))

solver.plot("Varying beta")
solver2.plot("Constant beta")

"""
run SIR_class.py
With beta varying between 0.0005 and 0.0001, the maximum number of infected is 745.41.

With beta at a constant 0.0005, the maximum number of infected is 897.871

-Output is two plots-
"""
