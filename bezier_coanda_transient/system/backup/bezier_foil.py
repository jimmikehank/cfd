# Code developed by James Henry for Bezier curve matching to airfoil coordinates.

# This function generates Bernstein polynomials of order 'n'
def bernstein(n,i,npoints):
    import numpy as np
    import math
    mappts = np.linspace(0,1,npoints)
    preterm = math.factorial(n)/(math.factorial(i)*math.factorial(n-i))
    poly = preterm * mappts**i * (1-mappts)**(n-i)
    bspoly = np.vstack([poly]).T
    return bspoly

# This function generates Bezier curves from the Bernstein polynomials
def bezier_curve(control_points,Nt):
    import numpy as np
    n = np.shape(control_points)[0]
    bezier = np.zeros([Nt,2])
    for i in range(n):
        bernpoly = bernstein(n-1,i,Nt)
        add = (bernpoly*control_points[i,:])
        bezier += add
    return bezier

    bezfoil = np.vstack([upper[::-1],lower])
    return bezfoil

# This function creates an airfoil from the Bezier curve, the coordinates
# follow a path from TE over top surface, to lower surface back to TE
def init_bezfoil(m, control_points_upper, control_points_lower = [],symmetric = True):
    import numpy as np
    inverse = np.array([1,-1])
    upper = bezier_curve(control_points_upper,m)
    if symmetric == True:
        lower = (upper * inverse)
    else:
        lower = bezier_curve(control_points_lower,m)
    foil = np.vstack([upper[::-1],lower])
    return foil

# This functions loads in the coordinates of the airfoil from the csv file
def load_airfoil(filename):
    import numpy as np
    coords = np.loadtxt(filename,delimiter=',')
    coords = (coords / np.max(coords[:,0]))
    return coords

# This function interpolates the loaded coordinates to match the x coordinates
# of the Bezier curves for error calculations
def interp_foil(bezfoil, airfoil):
    from scipy import interpolate
    import numpy as np
    bez_mids = np.where(bezfoil[:,0] == 0)[0]
    af_mid = np.where(airfoil[:,0] == 0)[0][0]
    airfoil_upper = airfoil[0:af_mid+1,:]
    airfoil_lower = airfoil[af_mid:,:]
    interp_upper_foil = interpolate.interp1d(airfoil_upper[:,0],airfoil_upper[:,1],kind='linear',fill_value='extrapolate')
    interp_lower_foil = interpolate.interp1d(airfoil_lower[:,0],airfoil_lower[:,1],kind='linear',fill_value='extrapolate')
    afint_upper = interp_upper_foil(bezfoil[0:bez_mids[1],0])
    afint_lower = interp_lower_foil(bezfoil[bez_mids[1]:,0])
    afint = np.hstack([afint_upper,afint_lower])
    return afint

# This function calculates the squared error of the airfoil
def bezier_error(bezfoil, airfoil):
    import numpy as np
#     error = 1 / len(bezfoil) * np.sum(np.abs((bezfoil[:,1] - airfoil)))
    error = 1 / len(bezfoil) * np.sum(((bezfoil[:,1] - airfoil)**2))
    return error

# This function calculates the gradient (de/dx) where e is the error and x is the control points
def find_grad(airfoil,control_points,m,control_points_lower=[],symmetric=True):
    import numpy as np
    dx = 1e-6
    err_grad = np.zeros([np.shape(control_points)[0],2])
    err_grad_lower = np.zeros([np.shape(control_points)[0],2])
    cp = control_points
    cpl = control_points_lower
    err_init = bezier_error(init_bezfoil(m,cp,cpl,symmetric),airfoil)
    if symmetric:
        for i in range(1,np.shape(control_points)[0]-1):
            if i == 1:
                pass
            else:
                cp = control_points
                cp[i,0] += dx
                newfoil = init_bezfoil(m,cp,cpl,symmetric)
                err_grad[i,0] = (bezier_error(newfoil,airfoil)-err_init)/dx
            cp = control_points
            cp[i,1] += dx

            newfoil = init_bezfoil(m,cp,cpl,symmetric)
            err_grad[i,1] = (bezier_error(newfoil,airfoil)-err_init)/dx
        err_grad_lower = err_grad*np.array([1,-1])
    else:
        for i in range(1,np.shape(control_points)[0]-1):
            if i == 1:
                pass
            else:
                cp = control_points
                cp[i,0] += dx
                newfoil = init_bezfoil(m,cp,cpl,symmetric)
                err_grad[i,0] = (bezier_error(newfoil,airfoil)-err_init)/dx
            cp = control_points
            cp[i,1] += dx
            newfoil = init_bezfoil(m,cp,cpl,symmetric)
            err_grad[i,1] = (bezier_error(newfoil,airfoil)-err_init)/dx
        for i in range(1,np.shape(control_points_lower)[0]-1):
            if i == 1:
                pass
            else:
                cpl = control_points_lower
                cpl[i,0] += dx
                newfoil = init_bezfoil(m,control_points,cpl,symmetric)
                err_grad_lower[i,0] = (bezier_error(newfoil,airfoil)-err_init)/dx
            cp = control_points_lower
            cp[i,1] += dx
            newfoil = init_bezfoil(m,control_points,cpl,symmetric)
            err_grad_lower[i,1] = (bezier_error(newfoil,airfoil)-err_init)/dx
    return err_grad, err_grad_lower

# This function uses all of the above code to optimize the bezier shape to that of the airfoil
def foil_opt(control_points,af_filename,chord = 1., eps=1e-6, deps=1e-12, m=101, step=1, debug=False,control_points_lower=[],sym=True):
    import numpy as np
    from matplotlib import pyplot as plt
    control_points = control_points
    airfoil_raw = load_airfoil(af_filename)
    bezfoil = init_bezfoil(m, control_points)
    airfoil = interp_foil(bezfoil,airfoil_raw)
    err = 1
    err_prev = .5
    derr = 1
    err_store = []
    iters = 0
    while err > eps and derr > deps:
        err = bezier_error(bezfoil,airfoil)
        err_store = np.append(err_store,err)
        grad, grad_lower = find_grad(airfoil,control_points,m,control_points_lower=control_points_lower,symmetric=sym)
        control_points -= grad*(step)
        if sym:
            pass
        else:
            control_points_lower -= grad_lower*step
        bezfoil = init_bezfoil(m, control_points,control_points_lower=control_points_lower,symmetric=sym)
        derr = np.abs(err_prev - err)
        err_prev = err
        iters+=1
        if iters % 100 == 0 and debug:
            print(err)
    if debug:
        plt.figure(figsize=[20,3])
        plt.plot(airfoil_raw[:,0],airfoil_raw[:,1],'b-')
        plt.plot(bezfoil[:,0],bezfoil[:,1],'r--')
        plt.axis('equal')
        plt.ylim([-.08,.08])
        plt.tick_params(axis='both',labelsize=16)
        plt.legend(['Original Airfoil','Optimized Airfoil'],fontsize=16)
        plt.savefig('airfoil_comparison.png')
    bezfoil = (bezfoil - control_points[-1,:]) * chord
    control_points = (control_points - control_points[-1,:]) * chord
    control_points_lower = (control_points_lower - control_points_lower[-1,:]) * chord
    if sym:
        control_points_lower = control_points * np.array([1,-1])
        return bezfoil, control_points, control_points_lower,iters
    else:
        return bezfoil, control_points, control_points_lower, iters


def coanda_foil(r,t,h,cpu,cpl=[],m=101):
    import numpy as np

    te_height = r+t+h
    foil_upper_init = bezier_curve(cpu,m)

    if cpl == []:
        cpl = cpu*np.array([1,-1])
    else:
        pass
    foil_lower_init = bezier_curve(cpl,m)
    intersect_upper = 0
    intersect_lower = 0
    for i in range(1,np.shape(foil_upper_init)[0]-1):
        if foil_upper_init[i-1,1] > te_height and foil_upper_init[i,1] <= te_height:
            intersect_upper = i
            break

    for i in range(np.shape(foil_lower_init)[0]-1):
        if foil_lower_init[i-1,1] < -te_height and foil_lower_init[i,1] >= -te_height:
            intersect_lower = i
            break
    eps = 1e-6
    coanda_pts =np.array([[0,0],[-r,r],[-r+eps,r+h],[-r+eps*2,r+h+t],[foil_upper_init[intersect_upper-1,0],foil_upper_init[intersect_upper-1,1]],[-r,-r],[-r+eps,-r-h],[-r+eps*2,-r-h-t],[foil_lower_init[intersect_lower-1,0],foil_lower_init[intersect_lower-1,1]]])
    foil_upper = foil_upper_init[:intersect_upper,:]
    foil_lower = foil_lower_init[:intersect_lower,:]

    return foil_upper, foil_lower, coanda_pts
