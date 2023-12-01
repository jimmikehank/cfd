from processing import *
import numpy as np
import os
import time


cmu = 0.025
rho = 1.17
mu  = 1.82e-5
c   = 0.3
b   = 0.156
parallel = False



if parallel:
    run_command = "mpirun -np 6 rhoSimpleFoam -parallel"
    decompose = "decomposePar"
    reconstruct = "reconstructPar -latestTime"
else:
    run_command = "rhoSimpleFoam"



def initializer(folder = './'):
    import os
    import numpy as np
    target_dir = '/home/james/Documents/research/cfd/bezier_coanda_transient/1'
    files = os.listdir()
    cleaned = [int(x) for x in files if check_float(x)]
    print(cleaned)
    ic_file = np.max(cleaned)
    os.system('cp -r ./{} {}'.format(ic_file, target_dir))


initializer()
