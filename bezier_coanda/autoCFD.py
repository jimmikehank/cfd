from processing import *
import numpy as np
import os
import time


cmu = 0.0075
rho = 1.17
mu  = 1.82e-5
c   = 1.0
b   = 0.156

re_range = np.linspace(0.25e6,4e6,16)
# re_range = np.linspace(2.25e6,3.5e6,6)
continue_bool = input("Reynolds number sweep: {}\n Proceed (y/n): ".format(re_range)).lower()
if continue_bool == 'y':
    for re in re_range:
        single_run(cmu, cmu, re, rho, mu, c, b, .005*cmu, urf=0.4, parallel=False)

elif continue_bool == 'n':
    print("\nSimulation exited\n")
else:
    print("Input must be y or n, simulation exited")
