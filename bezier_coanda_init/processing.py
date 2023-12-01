def max_velocity(folder='./', debug=False):
    import os
    import shutil
    import numpy as np
    files = os.listdir(folder)
    cleaned = [float(x) for x in files if check_float(x)]
    filename = np.max(cleaned)
    if int(filename) == filename:
        filename = '{}/{}/U'.format(folder,int(filename))
    else:
        filename = '{}/{}/U'.format(folder,filename)
    with open(filename) as f:
        Utext = f.readlines()
    numlines = int(Utext[19])
    start = 21
    end = start+numlines-1
    Umax = 0
    for i in range(start,end):
        line = Utext[i][1:-2]
        vals = [float(x) for x in line.split()]
        Umag = np.sqrt(vals[0]**2 + vals[1]**2)
        if Umag > Umax:
            Umax = Umag
    return Umax

def cmu_calc(U, mdot, c, b, rho, folder = './'):
    ujet = max_velocity(folder)
    cmu = mdot * ujet / (0.5 * rho * U**2 * c * b)
    return cmu

def cmu_feedback(cmu_command, cmu_actual):
    e = (cmu_command - cmu_actual)
    return e

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

def check_force(folder):
    import os
    import numpy as np
    files = os.listdir(folder)
    cleaned = [int(x) for x in files if check_float(x)]
    if os.path.exists(folder+"/postProcessing/forces/{}".format(np.max(cleaned))):
        return False
    else:
        return True

def cmu_openloop(cmu, U, c, b, rho):
    import numpy as np
    Q = 0.5 * rho * U ** 2 * c * b
    mdot = np.sqrt(cmu * Q / (150 / 0.007573))
    return mdot

def re_convert(re, rho, c, mu = 1.82e-5):
    U = re * mu / (rho * c)
    return U

def run_fun(parallel = False):
    import os
    if parallel:
        decompose = "decomposePar"
        run = "mpirun -np 6 rhoSimpleFoam -parallel"
        reconstruct = "reconstructPar -latestTime"
        os.system(decompose)
        os.system(run)
        os.system(reconstruct)
    else:
        run = "rhoSimpleFoam"
        os.system(run)


def single_run(cmu_target, cmu_command, re_command, rho, mu, c, b, eps, urf, parallel=False):
    import os
    import time
    U = re_convert(re_command, rho, c, mu)
    print(U)
    time.sleep(3)
    mdot = cmu_openloop(cmu_command, U, c, b, rho)
    mesh_command = "python3 autoMesh.py --clean true --mdot {} --meanFlow {} --airfoil NACA0015".format(mdot, U)
    print("AutoMesh command: {}".format(mesh_command))
    blockmesh = "blockMesh"
    run_command = "rhoSimpleFoam"
    os.system(mesh_command)
    os.system(blockmesh)
    run_fun(parallel)
    cmu_actual = cmu_calc(U, mdot, c, b, rho)
    print("Cmu Measured: {}".format(cmu_actual))
    e = cmu_feedback(cmu_target, cmu_actual)
    kp = urf
    print("Cmu error: {}".format(e))

    if abs(e) < eps:
        print("Single run completed")
        print("\n\tFinal Values:\n")
        print("Vjet: {}\nCmu: {}\nCmu error: {}".format(max_velocity('./'), cmu_actual, e))
        save_command = "python3 autoMesh.py --store True --runName /c{}/cmu{}/re{}".format(c,cmu_target,re_command)
        os.system(save_command)
        return 0
    else:
        print("Single run did not meet Cmu criteria, running again with Cmu = {}".format(cmu_command + kp*e))
        print("For diagnostics, Vjet = {}".format(max_velocity('./')))
        time.sleep(5)
        single_run(cmu_target, cmu_command + kp*e, re_command, rho, mu, c, b, eps, urf, parallel)
