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

def retrieve_lift(folder,debug=False):
    import os
    import shutil
    import numpy as np

    files = os.listdir(folder)
    cleaned = [x for x in files if check_float(x)]
    files = sorted(cleaned)
    needs_forces = check_force(folder)
    if debug:
        print(folder,needs_forces)
    if needs_forces:
        force_command = "rhoSimpleFoam -postProcess -case {} -func forces".format(folder)
        os.system(force_command)
    else:
        pass

    forces = np.zeros(3)
    moments = np.zeros(3)
    time = np.array([])

    for file in files:
        with open("{}/postProcessing/forces/{}/forces.dat".format(folder,file)) as f:
            full = f.readlines()
            line = full[3]
            starts = []
            ends = []
            for i in range(len(line)):
                if line[i] == '(' and line[i+1] != '(':
                    starts.append(i+1)
                elif line[i] == ')' and line[i-1] != ')':
                    ends.append(i)
            pressure_forces = np.array([float(x) for x in line[starts[0]:ends[0]].split()])
            viscous_forces = np.array([float(x) for x in line[starts[1]:ends[1]].split()])

            pressure_moments = np.array([float(x) for x in line[starts[2]:ends[2]].split()])
            viscous_moments = np.array([float(x) for x in line[starts[3]:ends[3]].split()])
            forces = np.vstack([forces, pressure_forces + viscous_forces])
            moments = np.vstack([moments, pressure_moments + viscous_moments])
            time = np.append(time,float(file))
    return forces, moments, time

def cmu_openloop(cmu, U, c, b, rho):
    import numpy as np
    Q = 0.5 * rho * U ** 2 * c * b
    mdot = np.sqrt(cmu * Q / (190 / 0.007573))
    return mdot

def re_convert(re, rho, c, mu = 1.82e-5):
    U = re * mu / (rho * c)
    return U

def run_fun(parallel = False):
    import os
    if parallel:
        decompose = "decomposePar"
        run = "mpirun -np 2 rhoSimpleFoam -parallel"
        reconstruct = "reconstructPar -latestTime"
        os.system(decompose)
        os.system(run)
        os.system(reconstruct)
    else:
        run = "rhoSimpleFoam"
        os.system(run)


def single_run(cmu_target, cmu_command, re_command, rho, mu, c, b, eps, urf, airfoil='naca0015',parallel=False):
    import os
    import time
    U = re_convert(re_command, rho, c, mu)
    print(U)
    time.sleep(3)
    mdot = cmu_openloop(cmu_command, U, c, b, rho)
    mesh_command = "python3 autoMesh.py --clean true --mdot {} --meanFlow {} --airfoil ".format(mdot, U, airfoil)
    print("AutoMesh command: {}".format(mesh_command))
    blockmesh = "blockMesh"
    run_command = "rhoSimpleFoam"
    os.system(mesh_command)
    os.system(blockmesh)
    run_fun(parallel)
    cmu_actual = cmu_calc(U, mdot, c, b, rho)
    print("Cmu Measured: {}".format(cmu_actual))
    e = cmu_feedback(cmu_target, cmu_actual)
    if cmu_target == 0:
        e = 0
    kp = urf
    print("Cmu error: {}".format(e))

    if abs(e) <= eps:
        print("Single run completed")
        print("\n\tFinal Values:\n")
        print("Vjet: {}\nCmu: {}\nCmu error: {}".format(max_velocity('./'), cmu_actual, e))
        save_command = "python3 autoMesh.py --store True --runName /c{}/cmu{}/re{}".format(c,round(cmu_target,4),re_command)
        # save_command = "python3 autoMesh.py --store True --runName /c{}/cmu{}".format(c,cmu_target,re_command)
        os.system(save_command)
        return 0
    else:
        print("Single run did not meet Cmu criteria, running again with Cmu = {}".format(cmu_command + kp*e))
        print("For diagnostics, Vjet = {}".format(max_velocity('./')))
        time.sleep(5)
        single_run(cmu_target, cmu_command + kp*e, re_command, rho, mu, c, b, eps, urf, parallel)
