"""Implementing a full-scale war on the zombies"""

import numpy as np
import SIZR
import ODESolver

def omega(t):
    """Defining omega"""
    s = 0
    T = [5, 10, 18]  # The tree t-values in which to activate omega

    if t >= 5 and t <= 7:
        a = 50*alpha
    elif t >= 10 and t <= 12:
        a = 50*alpha
    elif t >= 18 and t <= 20:
        a = 50*alpha
    else:
        a = 0

    for i in range(len(T)): # Omega is calculated using a loop
        s += np.exp(-0.5*((t-T[i])/lc_sig)**2)
    return a*s              # And returned with the a-value

class SIZR_war(SIZR.ProblemSIZR):
    def __init__(self, sigma, beta, dS, dI, rho, alpha, S0, I0, Z0, R0, T, omega):
        SIZR.ProblemSIZR.__init__(self, sigma, beta, dS, dI, rho, alpha, S0, I0, Z0, R0, T)
        self.omega = omega


        if isinstance(omega, (float,int)):
            """Allowing omega to be a function"""
            self.omega = lambda t: omega
        elif callable(omega):
            self.omega = omega


    def __call__(self, u, t):
        S, I, Z, R = u
        return [self.sigma(t) - self.beta(t)*S*Z - self.dS(t)*S,                 # S'(t)
                self.beta(t)*S*Z - self.rho(t)*I - self.dI(t)*I,                 # I'(t)
                self.rho(t)*I - (self.alpha(t)+self.omega(t))*S*Z,               # Z'(t) with omega(t)
                self.dS(t)*S + self.dI(t)*I + (self.alpha(t)+self.omega(t))*S*Z] # R'(t) with omega(t)

class SolverSIZR_war(SIZR.SolverSIZR):
    def __init__(self, problem, dt):
        SIZR.SolverSIZR.__init__(self, problem, dt)

    def solve(self, method=ODESolver.RungeKutta4):
            self.solver = method(self.problem)
            ic = [self.problem.S0, self.problem.I0, self.problem.Z0, self.problem.R0]
            self.solver.set_initial_condition(ic)
            n = int(round(self.problem.T/self.dt))
            t = np.linspace(0, self.problem.T, n+1)
            u, self.t = self.solver.solve(t)
            self.S, self.I, self.Z, self.R = u[:,0], u[:,1], u[:,2], u[:,3]
            return u, t

    def plot(self):
        SIZR.SolverSIZR.plot(self, "War on Zombies")



S0 = 50; Z0 = 3; I0 = 0; R0 = 0    # Initial conditions
beta = 0.03
alpha = 0.2*beta
lc_sig = 0.5                       # Lower-case sigma for calculating omega
dS = 0; dI = 0; sigma = 0; rho = 1
T = 20
problem = SIZR_war(sigma=sigma, beta=beta, rho=rho,
                               alpha=alpha, dI=dI, dS=dS, S0=S0,
                               I0=I0, Z0=Z0, R0=R0, T=T, omega=omega)

solver = SolverSIZR_war(problem, dt=0.25)
solver.plot()

"""We see at the first attack at t=5, the a large amount of zombies are killed,
however, there are also alot of human casulties. Therefore, in the next attack
at t = 10, there are fewer humans, and therefore not alot of zombies are killed.
The last attack at t = 18 is barely noticiable, due to the few remaining humans.
However, the amount of zombies has stabilized (for now), so maybe there is still hope"""

"""
run war_on_zombies.py
-Output is a plot-
"""
