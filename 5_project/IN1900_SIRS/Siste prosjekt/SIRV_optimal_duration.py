# Let the vaccination campaign last for Vt days
"""Compute the maximum number of infected people, maxt I.t/, as a
function of Vt, by running the model for
VT D 0; 1; 2 : : : ; 31. Plot this function.
Determine from the plot the optimal VT , i.e., the smallest
vaccination period VT such that increasing VT has negligible
effect on the maximum number of infected people."""

import numpy as np
import SIRV_varying_p as SIRV

Vt = np.linspace(0, 31, 31) # Vt-values
Imax = np.zeros(len(Vt))    # Empty array for maxmium amount of infected

for i in range(len(Vt)): # A loop that creates a new instance for each value of Vt

    problem = SIRV.SIRV_varying_p(beta=lambda t: 0.0005 if t <= 12 else 0.0001,
                               nu = 0.1, S0 = 1500, I0=1, R0=0, T=60, V0 = 0,
                               p=lambda t: 0.1 if t >= 6 and t <= 6+i else 0)

    solver = SIRV.SolverSIRV_varying_p(problem, dt=0.5)
    u = solver.solve()[0]                                # We extract the u-values
    Imax[i] = np.amax(u[:,1])                            # And add the max-values of the infected-column


import matplotlib.pyplot as plt
plt.plot(Vt, Imax, label = "Maximum number of infected")
plt.xlabel("Days with vaccination", fontsize = 20); plt.ylabel("Maximum infected", fontsize = 20)
plt.legend(fontsize=30)
plt.show()

"""We can see from the plot that starting with Vt = 6, the vaccination
has no negligible effect effect on the maximum number of infected people
as the curve is constant for Vt >= 6. This means the ideal t-value would be
6+Vt_ideal = 12"""

"""
run SIRV_optimal_duration.py
-Output is a plot-
"""
