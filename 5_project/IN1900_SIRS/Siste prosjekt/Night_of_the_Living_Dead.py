"""Simulating the three phases of Night of The Living Dead"""

import SIZR

t1 = 4     # The first 4 hours
t2 = 24+t1 # The next 24 hours
t3 = t2+5  # The last 5 hours

sigma  = lambda t: 20 if t <= t1 else 2 if t <= t2 and t > t1 else 0
beta   = lambda t: 0.03 if t <= t1 else 0.0012 if t <= t2 and t > t1 else 0
alpha  = lambda t: 0 if t <= t1 else 0.0016 if t <= t2 and t > t1 else 0.006
dI     = lambda t: 0.014 if t > t1 and t <= t2 else 0
dS     = lambda t: 0.0067 if t <= t3 and t > t2 else 0
rho = 1                          # Rho is constant throughout the time interval
S0 = 60; Z0 = 1; I0 = 0; R0 = 0  # Initial conditions
T = t3                           # Final t-value

"""We use the imported functions from the SIZR-module to solve the problem"""

problem = SIZR.ProblemSIZR(sigma=sigma, beta=beta, rho=rho, alpha=alpha,
                           dI=dI, dS=dS, S0=S0, I0=I0, Z0=Z0, R0=R0, T=T)
solver = SIZR.SolverSIZR(problem, dt=0.25)
solver.plot("The three phases of Night of The Living Dead")

"""
run Night_of_the_Living_Dead.py
-Output is a plot-
"""
