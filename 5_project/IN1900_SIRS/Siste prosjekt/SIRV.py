import SIR_class
import ODESolver
import numpy as np

class ProblemSIRV(SIR_class.ProblemSIR):
    def __init__(self, nu, beta, S0, I0, R0, T, V0, p):
        SIR_class.ProblemSIR.__init__(self, nu, beta, S0, I0, R0, T)
        self.V0, self.p = V0, p

    def __call__(self, u, t):
        """ODE's to solve"""
        S, I, R, V = u
        return [-self.beta(t)*S*I - self.p*S,   # New S'(t)
                self.beta(t)*S*I-self.nu(t)*I,  # I'(t)
                self.nu(t)*I,                   # R'(t)
                self.p*S]                       # V'(t) Vaccinated people


class SolverSIRV(SIR_class.SolverSIR):
    def __init__(self, problem, dt):
        SIR_class.SolverSIR.__init__(self, problem, dt)

    def solve(self, method=ODESolver.RungeKutta4):
        self.solver = method(self.problem)
        ic = [self.problem.S0, self.problem.I0, self.problem.R0, self.problem.V0] # A new initial condition V0 is added
        self.solver.set_initial_condition(ic)
        n = int(round(self.problem.T/self.dt))
        t = np.linspace(0, self.problem.T, n+1)
        u, self.t = self.solver.solve(t)
        self.S, self.I, self.R = u[:,0], u[:,1], u[:,2]
        return u, t

    def plot(self, title):
        """Plotting SIR with the addition of V"""
        import matplotlib.pyplot as plt
        u, t = self.solve()
        V = u[:,3]
        plt.plot(t, V, label = "Vaccinated")
        SIR_class.SolverSIR.plot(self, title)

problem = ProblemSIRV(beta=lambda t: 0.0005 if t <= 12 else 0.0001,
                      nu = 0.1, S0 = 1500, I0=1, R0=0, T=60, V0 = 0, p= 0.1)
solver = SolverSIRV(problem, 0.5)
solver.plot("SIR with vaccination")


"""We see that as the rate of vaccinated people increases, the amount
of susceptibles decreases. Not because they get infected (we see this as the
infected-curve never reaches a high point), but because when they are vaccinated,
they are no longer susceptible to the disease."""

"""
run SIRV.py
-Output is a plot-
"""
