import numpy as np
import os
from processing import *
from matplotlib import pyplot as plt
import time

meanflow = 30
S = 0.156 * 0.3
c = 0.3
rho = 1.17
QS = 0.5*1.17*meanflow**2*S
dP = 5000
press_base = 100000 + dP/2
N = 41
dPs = np.linspace(-dP/2,dP/2,N)
file = './output/p_vs_f.txt'

if os.path.exists(file):
    pass
else:
    np.savetxt(file,np.array([]))

def runone(mesh,save='empty'):
    import os
    block = 'blockMesh'
    run = 'rhoSimpleFoam'

    os.system(mesh)
    os.system(block)
    os.system(run)
    # os.system(save)
    forces, moments, _time = retrieve_lift('./')
    return forces

re = 625000
meanflow = re_convert(re,rho,c)
aoas = np.arange(0.0,10.1,2.5)
mdots = np.arange(0,0.0061,0.0005)
i = 0
for aoa in aoas:
    for mdot in mdots:
        pressure = mdot_to_pressure_bc(mdot)
        mesh = 'python3 autoMesh.py --clean True --airfoil naca0015 --meanFlow {} --pressUp {} --aoa {}'.format(meanflow,pressure,aoa)
        block = 'blockMesh'
        run = 'rhoSimpleFoam'
        os.system(mesh)
        os.system(block)
        os.system(run)
        mdot_actual = massflow_bc('./')
        save = 'python3 autoMesh.py --store True --runName aoa{}_mdot{}'.format(aoa,mdot)
        print(save)
        os.system(save)
