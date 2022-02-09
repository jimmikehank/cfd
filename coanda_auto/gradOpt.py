import numpy as np
import os

xinit = np.array([
    0.003,    # Rc - Coanda cylinder radius
    0.00015,  # te - Slot exit height
    0.00018,  # tu - Upper surface thickness at exit
    0.5,      # ru - Contour radius of upper surface to exit
    60.0,     # ai - internal contraction angle
])

x = xinit
mdot = 0.05

def meshCommand(x):
    command = 'python3 autoMesh.py -Rc {} -te {} -tu {} -ru {} -ai {}'.format(x[0],x[1],x[2],x[3],x[4])
    return command

def bocon(mdot):
    command = 'python3 autoBCu.py {}'.format(mdot)
    return command

def cleanhouse(iter):
    retain = ['0', 'constant', 'system','autoCFD.py','converged','test_bench','autoBC.py','autoBCu.py','autoMesh.py','images','gradOpt.py']
    dirs = os.listdir()

    for item in dirs:
        if item not in retain:
            if check_float(item):
                if float(item)%1 == 0:
                    it = int(item)
                else:
                    it = float(item)
                boop = boop + [it]
            else:
                shutil.rmtree(item)

    new_converged = str(max(boop))

    foldname = 'iteration_{}'.format(iter)
    os.system('mkdir converged/{}'.format(foldname))
    os.system('mv {} converged/{}'.format(new_converged,foldname))

    for item in copy:
        os.system('cp {} converged/{} -r'.format(item,foldname))

    dirs = os.listdir()

    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)

    return foldname

def retrieveLift():
    import os
    os.system('rhoSimpleFoam -postProcess -func forces -latestTime')
    forceFold = os.listdir('./postProcessing/forces')
    forceFold = forceFold[0]
    sel = []
    with open('./postProcessing/forces/{}/forces.dat'.format(forceFold)) as f:
        data = f.readlines()[3]
    for i in range(len(data)):
        if data[i] == '(' and data[i+1] == '(':
            j = i
            g = 0
            while g < 2:
                if data[j] == ' ':
                    sel.append(j)
                    g = g + 1
                    j = j+1
                else:
                    j = j+1

            break
    lift = float(data[sel[0]:sel[1]])
    return lift

def singlerun(x,mdot):
    import os
    commandMesh = meshCommand(x)
    commandBC = bocon(mdot)

    os.system(commandMesh)
    os.system(commandBC)
    os.system('decomposePar')
    os.system('mpirun -np 4 rhoSimpleFoam -parallel')
    os.system('reconstructPar -latestTime')
