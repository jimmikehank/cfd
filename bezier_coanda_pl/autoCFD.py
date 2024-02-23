import numpy as np
import os
from processing import *
from matplotlib import pyplot as plt
import time

meanflow = 30
S = 0.156 * 0.3
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
#
#
# for i in range(N):
#     Pu = press_base + dPs[i]
#     Pl = press_base - dPs[i]
#     if i == 0:
#         print('Running with Pu: {} and Pl : {}'.format(Pu,Pl))
#         input("Continue?")
#
#     mesh_naca = 'python3 autoMesh.py --clean True --airfoil naca0015 --meanFlow {} --pressUp {} --pressLo {}'.format(meanflow,Pu,Pl)
#     save = 'python3 autoMesh.py --store True --runName pu{}_pl{}'.format(Pu,Pl)
# # mesh_naca0 = 'python3 autoMesh.py --clean True --airfoil naca0015 --meanFlow {} --pressUp {}'.format(meanflow,100000)
# # mesh_opt1 = 'python3 autoMesh.py --clean True --airfoil opt1 --meanFlow {} --pressUp {}'.format(meanflow,pressUp)
# # mesh_opt0 = 'python3 autoMesh.py --clean True --airfoil opt1 --meanFlow {} --pressUp {}'.format(meanflow,100000)
#
#     forces = runone(mesh_naca,save)
#     drag = forces[-1,0] / (0.5 * 1.17 * meanflow ** 2 * S)
#     lift = forces[-1,1] / (0.5 * 1.17 * meanflow ** 2 * S)
#
#     os.system('foamToVTK -latestTime')
#     U = max_velocity('./')
#     mdot = massflow_bc('./')
#     T = 2 * U * mdot
#     point = np.array([Pu,Pl,lift,drag])
#     if i == 0:
#         np.savetxt(file,point)
#     else:
#         data = np.loadtxt(file)
#         data = np.vstack([data,point])
#         np.savetxt(file,data)

# print("Max velocity: {}\nThrust AFC: {}\nMass flow rate: {}\nRaw Drag: {}\nDrag plus thrust: {}".format(U, T, mdot, drag, drag+T))


mesh_naca0 = 'python3 autoMesh.py --clean True --airfoil naca0015 --meanFlow {} --pressUp {}'.format(meanflow,115000)
forces = runone(mesh_naca0)

drag = forces[-1,0] / (0.5 * 1.17 * meanflow ** 2 * S)
lift = forces[-1,1] / (0.5 * 1.17 * meanflow ** 2 * S)
os.system('foamToVTK')
mdot = massflow_bc('./')

print(lift,drag,mdot)

# forces = runone(mesh_naca)
# lift_naca = forces[-1,1]
# drag_naca = forces[-1,0]
#
# forces = runone(mesh_naca0)
# lift_naca0 = forces[-1,1]
# drag_naca0 = forces[-1,0]
#
# forces = runone(mesh_opt1)
# lift_opt1 = forces[-1,1]
# drag_opt1 = forces[-1,0]
#
# forces = runone(mesh_opt0)
# lift_off = forces[-1,1]
# drag_off = forces[-1,0]
#
# dragNACA = drag_naca - drag_naca0
# dragOPT = drag_opt1 - drag_off
# print(dragNACA,dragOPT)





# mesh = 'python3 autoMesh.py --clean true --airfoil naca0015 --meanFlow 30'
# block = "blockMesh"
# decomp = "decomposePar"
# runpar = "mpirun -np 4 rhoSimpleFoam -parallel"
# recon = "reconstructPar"
# run = "rhoSimpleFoam"
# force = "rhoSimpleFoam -postProcess -func forces"
# VTK = "foamToVTK -latestTime"
#
# def plenum_test():
#     mesh = 'python3 autoMesh.py --clean true --airfoil naca0015 --meanFlow 30 --iter 5'
#     block = "blockMesh"
#     run = "rhoSimpleFoam"
#     force = "rhoSimpleFoam -postProcess -func forces"
#
#     os.system(mesh)
#     os.system(block)
#     os.system(run)
#     os.system("foamToVTK -latestTime")
#     os.system(force)
#     f, m, tt = retrieve_lift('./')
#     L = f[-1,1]
#     mdot = massflow_bc()
#     Ujet = max_velocity()
#     print(mdot*Ujet / (QS))
#     print(L/QS)
#
# # plenum_test()
#
# ###### Optimization of Control Effectors ######
#
# save_file = '/home/james/Documents/research/completed_cases/coanda_opt/'
# max_iters = 25
# delta = 2e-7
# gamma = np.array([5e-11,5e-12,5e-11])*0.1
#
# dL = np.zeros([max_iters,4])
# dD = np.zeros([max_iters,4])
#
# # Initial values of R, h, t
# r = 0.12*0.0254
# h = 0.009*0.0254
# t = 0.012*0.0254
#
#
# ind = np.arange(1,max_iters+1,1)
#
# for i in range(max_iters):
#     if i == 0:
#         x = np.array([r,h,t])
#     else:
#         x = np.loadtxt(save_file +'{}.txt'.format(i))
#     for j in range(4):
#         mesh = 'python3 autoMesh.py --clean true --iter {} --airfoil naca0015 --meanFlow {} --control {} --delta {}'.format(i, meanflow, j, delta)
#         os.system(mesh)
#         os.system(block)
#         os.system(decomp)
#         os.system(runpar)
#         os.system(recon)
#         os.system(force)
#         os.system(VTK)
#         Ujet = max_velocity()
#         mdot = massflow_bc()
#         forces, moments, tt = retrieve_lift('./')
#
#         L = forces[-1,1] / (Ujet * mdot)
#         D = forces[-1,0] / (Ujet * mdot)
#         if j == 0:
#             dL[i,j] = L
#         else:
#             dL[i,j] = (L - dL[i,0])/delta
#
#     print("Gradient dL/dx found: {}".format(dL[i,1:]))
#     xdot = x + gamma * dL[i,1:]
#     np.savetxt('{}{}.txt'.format(save_file,i+1),xdot)
#
# np.savetxt('{}/sens.txt'.format(save_file),dL[:,1:])
# print("Final Values:\nCmu: {}\nCL: {}".format(mdot * Ujet / QS, L * mdot * Ujet / QS))
#
# plt.figure(figsize=[10,10])
# for i in range(3):
#     plt.plot(ind,dL[:,1+i])
# plt.xlabel('Iteration')
# plt.ylabel('Cost Function')
# plt.legend(["R","h","t"])
# plt.savefig('{}/gradient_plot.png'.format(save_file))
