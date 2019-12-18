"""Let the vaccination campaign start 6 days after the outbreak
of the disease, and let it last for 10 days"""
import SIRV
import SIR_class as SIR
import ODESolver
import numpy as np

class SIRV_varying_p(SIRV.ProblemSIRV):
    def __init__(self, nu, beta, S0, I0, R0, T, V0, p):
        SIRV.ProblemSIRV.__init__(self, nu, beta, S0, I0, R0, T, V0, p)

        if isinstance(p, (float, int)):
            self.p = lambda t: p
        elif callable(p):
            self.p = p

    def __call__(self, u, t):
        S, I, R, V = u
        return [-self.beta(t)*S*I - self.p(t)*S,
                self.beta(t)*S*I-self.nu(t)*I,
                self.nu(t)*I,
                self.p(t)*S]

class SolverSIRV_varying_p(SIRV.SolverSIRV):
    def __init__(self, problem, dt):
        SIRV.SolverSIRV.__init__(self, problem, dt)

    def solve(self, method=ODESolver.RungeKutta4):
        self.solver = method(self.problem)
        ic = [self.problem.S0, self.problem.I0, self.problem.R0, self.problem.V0]
        self.solver.set_initial_condition(ic)
        n = int(round(self.problem.T/self.dt))
        t = np.linspace(0, self.problem.T, n+1)
        u, self.t = self.solver.solve(t)
        self.S, self.I, self.R = u[:,0], u[:,1], u[:,2]
        return u, t


    def plot(self):
        SIRV.SolverSIRV.plot(self, "SIRV with varying p")


dt = 0.5
problem = SIRV_varying_p(beta=lambda t: 0.0005 if t <= 12 else 0.0001,
                           nu = 0.1, S0 = 1500, I0=1, R0=0, T=60, V0 = 0,
                           p=lambda t: 0.1 if t >= 6 and t <= 15 else 0)

solver = SolverSIRV_varying_p(problem, dt)
solver.plot()

"""
run SIRV_varying_p.py
-Output is a plot-
"""
