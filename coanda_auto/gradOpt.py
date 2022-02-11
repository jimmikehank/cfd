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

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

def cleanhouse(iter,keep):
    import os
    import shutil

    retain = ['0', 'constant', 'system', 'autoCFD.py', 'converged', 'test_bench', 'autoBC.py', 'autoBCu.py', 'autoMesh.py', 'images', 'gradOpt.py','output']
    copy = ['0', 'constant', 'system']

    dirs = os.listdir()
    boop = []
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

    if keep == True:
        foldname = 'iteration_{}'.format(iter)
        os.system('mkdir converged/{}'.format(foldname))
        os.system('mv {} converged/{}'.format(new_converged,foldname))

        for item in copy:
            os.system('cp {} converged/{} -r'.format(item,foldname))

        dirs = os.listdir()

        for item in dirs:
            if item not in retain:
                shutil.rmtree(item)
    else:

        dirs = os.listdir()

        for item in dirs:
            if item not in retain:
                shutil.rmtree(item)

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
    os.system('blockMesh')
    os.system('decomposePar')
    os.system('mpirun -np 8 rhoSimpleFoam -parallel')
    os.system('reconstructPar -latestTime')
    lift = retrieveLift()
    return lift

def find_gradient(x,mdot,iteration,selection):
    import numpy as np
    # Start with flow solution around initial x case
    j = iteration
    if j == 0:
        x_init = x
        sgn = np.ones(5)
    else:
        x_init = x[j,:]
        sgn = np.sign(x[j,:]-x[j-1,:])
    print('{} find_gradient -> x_init'.format(x_init))
    lift_init = singlerun(x_init,mdot)
    cleanhouse(iteration,True)

    # Compute gradient by modifying individual x values and solving flow over and over and over again.
    dLdX = np.zeros(np.size(x_init))
    xnew = np.zeros(np.size(x_init))

    for i in range(np.size(x_init)):
        if selection[i] == 1:
            xnew[:] = x_init[:]
            xnew[i] = x_init[i] * (1 + 0.01 * sgn[i])
            dX = xnew[i] - x_init[i]
            lift_i = singlerun(xnew,mdot)
            cleanhouse(iteration,False)
            dL = lift_i - lift_init
            with open('./output/debug.txt','a') as f:
                f.writelines('dL\t{}\t\tdX:\t{}\n'.format(dL,dX))
            dLdX[i] = dL/dX
        else:
            dLdX[i] = 0
    return(dLdX)


def update_design(x,dLdX,iteration):
    # Step forward to x(i+1) using gradient "ascent"
    import numpy as np
    if iteration == 0:
        gamma = np.sign(dLdX)
        dx = 1.01*x - x
        x = np.vstack([x,x + gamma * dx])
    else:
        res_x = x[j,:] - x[j-1,:]
        res_gradL = dLdX[j,:] - dLdX[j-1,:]
        # For non zero cases apply Barzilai-Borwein method for gamma - guarantees convergence to a local minimum
        gamma = np.append(gamma,np.dot(res_x,res_gradL)/np.dot(res_gradL,res_gradL))
        x = np.vstack([x, x[iteration-1,:] + (gamma[iteration] * dLdX[iteration,:])])

    return x



eps = 1
convergence = 0.1
iteration = 0
selection = np.ones(5)
# Actual optimization loop, convergence criteria subject to change to increase or decrease sensitivity as necessary
while eps > convergence:
    if iteration == 0:
        dLdX = find_gradient(x,mdot,iteration,selection)
        x = update_design(x,dLdX,iteration)
        eps = np.max(np.abs(dLdX))
    else:
        dLdX = np.vstack([dLdX,find_gradient(x,mdot,iteration,selection)])
        x = update_design(x,dLdX,iteration)
        eps = np.max(dLdX[iteration,:])
        for i in range(np.size(dLdX[iteration,:])):
            if np.abs(dLdX[iteration,i]) < convergence:
                selection[i] = 0
            else:
                selection[i] = 1
    iteration = iteration + 1
    np.savetxt('./output/X.txt',x)
    np.savetxt('./output/dLdX.txt',dLdX)

print('Solution Converged!')
