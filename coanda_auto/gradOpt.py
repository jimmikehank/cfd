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

def inprogress(dirs):
    import numpy as np
    import os
    iternumbers = np.array([])
    for item in dirs:
        for j in range(len(item)):
            if item[j] == '_':
                iternumbers = np.append(iternumbers,int(item[j+1:]))
    current = int(np.max(iternumbers))
    return current

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
        export_fold = '/home/james/Documents/research/converged_cases'
        foldname = 'iteration_{}'.format(iter)
        os.system('mkdir {}/{}'.format(export_fold,foldname))
        os.system('mv {} {}/{}'.format(new_converged,export_fold,foldname))

        for item in copy:
            os.system('cp {} {}/{} -r'.format(item,export_fold,foldname))

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
            xnew[i] = x_init[i] * (1 + 0.0004 * sgn[i])
            dX = xnew[i] - x_init[i]
            lift_i = singlerun(xnew,mdot)
            cleanhouse(iteration,False)
            dL = lift_i - lift_init
            dLdX[i] = dL/dX
        else:
            dLdX[i] = 0
    return(dLdX)

# def delta_x_calc(x,dLdX,iteration,flip_scaler,lim,gamma):
#     xnew = x[iteration,:]*(1 + (gamma * dLdX[iteration,:] / x[0,:]))
#     xout = np.zeros(np.size(xnew))
#     for i in range(np.size(xnew)):
#         if np.abs(xnew[i]-x[iteration,i])/x[iteration,i] > lim:
#             xout[i] = (1 + np.sign(xnew[i]-x[iteration,i]) * lim) * x[iteration,i]
#             with open('./output/debug.txt','a') as f:
#                 f.writelines("{}".format(xout[i]))
#         else:
#             xout[i] = xnew[i]
#             with open('./output/debug.txt','a') as f:
#                 f.writelines("{}\t".format(xout[i]))
#     with open('./output/debug.txt','a') as f:
#         f.writelines("\n")
#     return(xout)

def adam(x,dLdX,m,v,iteration):
    # Step forward in x using Adam algorithm
    import numpy as np
    order = np.array([1e-2, 2e-3, 5e-3, 1, 1])
    alpha = 0.001
    beta1 = 0.9
    beta2 = 0.999
    eps = 1e-8
    m = np.vstack([m, beta1 * m[iteration-1,:] + (1 - beta1)*dLdX[iteration,:]])
    v = np.vstack([v, beta2 * v[iteration-1,:] + (1 - beta2)*dLdX[iteration,:]**2])
    mhat = m[iteration,:]/(1-beta1)
    vhat = v[iteration,:]/(1-beta2)
    x = np.vstack([x, x[iteration,:] + order * alpha * mhat / (np.sqrt(vhat + eps))])
    return x, m, v


def update_design(x,dLdX,iteration,gamma, m, v):
    # Step forward to x(i+1) using gradient "ascent"
    import numpy as np
    j = iteration
    if iteration == 0:
        gamma = np.sign(dLdX)
        dx = 1.0004*x - x
        x = np.vstack([x,x + gamma * dx])
        flip_scaler = np.ones(5)
    else:
        res_x = (x[j,:] - x[j-1,:])/np.min(x[j,:] - x[j-1,:])
        res_gradL = dLdX[j,:] - dLdX[j-1,:]
        flip_scaler = np.zeros(5)
        for i in range(5):
            if np.sign(dLdX[j,i]) != np.sign(dLdX[j-1,i]):
                flip_scaler[i] = 0.2
            else:
                flip_scaler[i] = 1


        # Use Adam to step forward in x (gradient step switched from negative to positive)
        x, m, v = adam(x,dLdX,m,v,iteration)

    return x, gamma, flip_scaler, m, v

# Initialization loop - Checks for content in converged folder, then continues from that iteration

checkdir = os.listdir('/home/james/Documents/research/converged_cases/')
print(checkdir)
if checkdir == []:
    selection = np.array([1,1,1,0,0])
    iteration = 0
    eps = 1
else:
    selection = np.ones(5)
    iteration = inprogress(checkdir)+1
    x = np.loadtxt('./output/X.txt')
    dLdX = np.loadtxt('./output/dLdX.txt')
    eps = np.max(np.abs(dLdX[iteration-1]))
    flip = np.loadtxt('./output/Flip.txt')
    m = np.loadtxt('./output/m.txt')
    v = np.loadtxt('./output/v.txt')
    for j in range(5):
        k = dLdX[int(iteration-1),j]
        if k == 0:
            selection[j] = 0
            
print('Starting at iteration {} with variables {}.\n'.format(iteration,selection))
yes = input('If this is incorrect, enter X:  ')
if yes.lower() == 'x':
    exit()
convergence = 0.1
max_epochs = 300
gamma = 1

# Actual optimization loop, convergence criteria subject to change to increase or decrease sensitivity as necessary

while eps > convergence and iteration < max_epochs:
    if iteration == 0:
        m = np.zeros([1,5])
        v = np.zeros([1,5])
        print(m)
        dLdX = find_gradient(x,mdot,iteration,selection)
        x, gamma, flip, m, v = update_design(x,dLdX,iteration, gamma, m, v)
        eps = np.max(np.abs(dLdX))
        gamma = np.array([0])
    else:
        dLdX = np.vstack([dLdX,find_gradient(x,mdot,iteration,selection)])
        x, gamma, flipi, m, v = update_design(x,dLdX,iteration, gamma, m, v)
        flip = np.vstack([flip,flipi])
        eps = np.max(np.abs(dLdX[iteration,:]))
        for i in range(np.size(dLdX[iteration,:])):
            if np.abs(dLdX[iteration,i]) < convergence:
                selection[i] = 0
                with open('./output/selection_checker.txt','a') as f:
                    f.writelines('Iteration {}\nVariable {} converged with partial derivative: {}\n\n'.format(iteration,i,dLdX[iteration,i]))
            else:
                selection[i] = 1
    iteration = iteration + 1
    np.savetxt('./output/Flip.txt',flip)
    np.savetxt('./output/X.txt',x)
    np.savetxt('./output/dLdX.txt',dLdX)
    np.savetxt('./output/m.txt',m)
    np.savetxt('./output/v.txt',v)

print('Solution Converged!')
