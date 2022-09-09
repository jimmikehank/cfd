def bernstein(n,i,npoints):
    import numpy as np
    import math
    mappts = np.linspace(0,1,npoints)
    preterm = math.factorial(n)/(math.factorial(i)*math.factorial(n-i))
    poly = preterm * mappts**i * (1-mappts)**(n-i)
    bspoly = np.vstack([poly]).T
    return bspoly

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


def init_bezfoil(m, control_points_upper, control_points_lower = [],symmetric = False):
    import numpy as np
    inverse = np.array([1,-1])
    upper = bezier_curve(control_points_upper,m)
    lower = (upper * inverse)
    foil = np.vstack([upper[::-1],lower])
    return foil

def load_airfoil(filename):
    import numpy as np
    coords = np.loadtxt(filename,delimiter=',')
    coords = (coords / np.max(coords[:,0]))
    return coords

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

def bezier_error(bezfoil, airfoil):
    import numpy as np
#     error = 1 / len(bezfoil) * np.sum(np.abs((bezfoil[:,1] - airfoil)))
    error = 1 / len(bezfoil) * np.sum(((bezfoil[:,1] - airfoil)**2))
    return error

def find_grad(airfoil,control_points,m):
    import numpy as np
    err_grad = np.zeros([np.shape(control_points)[0],2])
    cp = control_points
    err_init = bezier_error(init_bezfoil(m,cp),airfoil)
    for i in range(1,np.shape(control_points)[0]-1):
        if i == 1:
            pass
        else:
            cp = control_points
            cp[i,0] += 1e-6
            newfoil = init_bezfoil(m,cp)
            err_grad[i,0] = (bezier_error(newfoil,airfoil)-err_init)/1e-6
        cp = control_points
        cp[i,1] += 1e-6
        newfoil = init_bezfoil(m,cp)
        err_grad[i,1] = (bezier_error(newfoil,airfoil)-err_init)/1e-6
    return err_grad

def foil_opt(control_points,af_filename,chord = 1., eps=1e-6,m=101, debug = False):
    import numpy as np
    control_points = control_points
    airfoil_raw = load_airfoil(af_filename)
    bezfoil = init_bezfoil(m, control_points)
    airfoil = interp_foil(bezfoil,airfoil_raw)
    err = 1
    err_store = []
    iters = 0
    while err > eps:
        err = bezier_error(bezfoil,airfoil)
        err_store = np.append(err_store,err)
        grad = find_grad(airfoil,control_points,m)
        control_points -= grad*(1)
        bezfoil = init_bezfoil(m, control_points)
        iters+=1
        if iters % 100 == 0 and debug:
            print(err)
    bezfoil = (bezfoil - control_points[-1,:]) * chord
    control_points = (control_points - control_points[-1,:]) * chord

    return bezfoil, control_points, iters

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
