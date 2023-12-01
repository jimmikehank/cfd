# Code developed by James Henry for Bezier curve matching to airfoil coordinates.
import os
import numpy as np
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

def interp_surface(bezier, foil_surf):
    from scipy import interpolate
    import numpy as np
    interper = interpolate.interp1d(foil_surf[:,0], foil_surf[:,1], kind='linear',fill_value='extrapolate')
    foil_surf_int = interper(bezier[:,0])
    return foil_surf_int

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

def find_grad_surf(surf_int,control_points,m):
    import numpy as np
    dx = 1e-8
    err_grad = np.zeros([np.shape(control_points)[0],2])
    err_init = bezier_error(bezier_curve(control_points,m),surf_int)

    for i in range(1,np.shape(control_points)[0]-1):
        if i == 1:
            pass
        else:
            cp = control_points
            cp[i,0] += dx
            newsurf = bezier_curve(cp,m)
            err_grad[i,0] = (bezier_error(newsurf,surf_int)-err_init)/dx

        cp = control_points
        cp[i,1] += dx
        newsurf = bezier_curve(cp,m)
        err_grad[i,1] = (bezier_error(newsurf,surf_int)-err_init)/dx

    return err_grad

def foil_opt_vec(control_points, af_vec, chord = 1., eps=1e-6, deps=1e-12, m=101, step=1, debug=False, control_points_lower=[],sym=True,max_iters=1e6):
    import numpy as np
    from matplotlib import pyplot as plt
    control_points = control_points
    bezfoil = init_bezfoil(m, control_points,control_points_lower=control_points_lower,symmetric=sym)
    airfoil = interp_foil(bezfoil,af_vec)
    err = 1
    err_prev = .5
    derr = 1
    err_store = []
    iters = 0
    while err > eps and derr > deps and iters < max_iters:
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
    bezfoil = (bezfoil - control_points[-1,:]) * chord
    control_points = (control_points - control_points[-1,:]) * chord
    control_points_lower = (control_points_lower - control_points_lower[-1,:]) * chord
    if sym:
        control_points_lower = control_points * np.array([1,-1])
        return bezfoil, control_points, control_points_lower,iters
    else:
        return bezfoil, control_points, control_points_lower, iters

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

def surface_opt(control_points, surf, m = 101, eps = 1e-6, deps = 1e-12, max_iters = 1e6, step = 1, debug=False):
    import numpy as np
    bezier = bezier_curve(control_points, m)
    err = 1
    err_prev = .5
    derr = 1
    err_store = []
    iters = 0
    while err > eps and derr > deps and iters < max_iters:
        surf_int = interp_surface(bezier,surf)
        err = bezier_error(bezier, surf_int)
        err_store = np.append(err_store,err)
        grad = find_grad_surf(surf_int, control_points, m)
        control_points -= grad * step
        iters+=1
        if iters % 1000 == 0 and debug:
            print(err)
        bezier = bezier_curve(control_points, m)
        derr = np.abs(err - err_prev)
        err_prev = err
    return bezier, control_points, iters

def read_sens_map(filename):
    import numpy as np
    with open(filename, 'r') as f:
        sens_raw = f.readlines()
        sens_raw = sens_raw[27:-4]
        N = len(sens_raw)
        sens_map = np.zeros([N,3])
        for i in range(N):
            line = sens_raw[i].replace('(','').replace(')','')
            line_float = np.array([])
            for item in line.split(' '):
                line_float = np.append(line_float, float(item))
            sens_map[i,:] = line_float
    return sens_map

def get_updated_shape(surface_name,sensitivity_map,scale=1,casedir=os.getcwd(),num_objectives=1):
    import pyvista as vtki
    import numpy as np
    target = casedir + '/VTK/bezier_coanda_opt_0/boundary/{}.vtp'.format(surface_name)
    print(target)
    wing = vtki.PolyData(target)
    pts_raw = wing.points
    M = np.shape(pts_raw)[0]
    N = int(M/2)
    t = np.linspace(0,4*np.pi,N)
    window = (1 - np.cos(t))/2
    window = np.vstack([window,window]).T
    pts = np.zeros([N,2])
    j = 0

    for i in range(M):
        if pts_raw[i,2] > 0:
            pts[j,:] = pts_raw[i,:2]
            j = j+1
        else:
            pass

    half_point = int(N/2)

    if num_objectives > 1:
        newpts = pts
        for k in range(num_objectives):
            try:
                newpts = newpts + sensitivity_map[:,:2,k]*scale[k]*window
            except TypeError:
                print("Either sensitivity map or scaling is not same dimension as number of objectives")
                break
    else:
        newpts = pts + sensitivity_map[:,:2]*scale*window
    afU = newpts[:half_point+1,:]
    afL = newpts[half_point+1:,:]
    afU = afU[::-1]

    airfoil_update = np.vstack([afU,afL]) - [np.min(newpts[:,0]),0]
    return airfoil_update

def get_mesh_data(vtk_target,sens):
    import pyvista as vtki
    raw_data = vtki.PolyData(vtk_target)
    faces_raw = raw_data.faces
    pts_raw = raw_data.points
    N = np.shape(sens)[0]
    M = np.shape(pts_raw)[0]
    faces = np.reshape(faces_raw,[N,5])
    pts = np.zeros([int(M/2),2])
    j = 0
    for i in range(M):
        if pts_raw[i,2] > 0:
            pts[j,:] = pts_raw[i,:2]
            j = j+1
        else:
            pass
    face_points = pts_raw[faces[:,1:3]]
    return pts, faces, face_points


def window_results(points, points_mod, blocks_x_L, scale = 1):
    import numpy as np
    corner_index = blocks_x_L * 2
    N = np.shape(points)[0]
    window = np.zeros([N,2])
    t = np.linspace(0,4*np.pi,corner_index)
    window[:corner_index,0] = (1-np.cos(t))/2
    window[:corner_index,1] = (1-np.cos(t))/2

    new_pts = points + window * points_mod * scale
    upper = new_pts[:blocks_x_L,:]
    lower = new_pts[blocks_x_L+1:2*blocks_x_L+1,:]
    lower = lower[::-1,:]
    return upper, lower

def coanda_surfs(r,t,h,cpu,cpl=[],m=101):
    import numpy as np

    te_height = r+t+h
    foil_upper_init = bezier_curve(cpu,m)

    if len(cpl) == []:
        cpl = cpu*np.array([1,-1])
    else:
        pass

    foil_lower_init = bezier_curve(cpl,m)
    eps = 1e-6
    coanda_pts =np.array([[0,0],[-r,r],[-r+eps,r+h],[-r+eps*2,r+h+t],foil_upper_init[-1,:],[-r,-r],[-r+eps,-r-h],[-r+eps*2,-r-h-t],[foil_lower_init[-1,0],foil_lower_init[-1,1]]])
    return coanda_pts



def coanda_foil(r,t,h,cpu,cpl=[],m=101):
    import numpy as np

    te_height = r+t+h
    foil_upper_init = bezier_curve(cpu,m)

    if len(cpl) == 0:
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

def find_gaps(wall_points):
    try:
        assert(np.shape(wall_points)[1] == 2)
    except AssertionError:
        print('Walls must be pairs of 2 points')
        return 0
    gaps = np.array([])
    for k in range(np.shape(wall_points)[0]):
        if wall_points[k,0].all() != wall_points[k-1,1].all():
            gaps = np.append(gaps,k)
    return gaps


def length(point1, point2):
    import numpy as np
    length = np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
    return length

def scale_points(faces,face_points,sens_map_wall):
    import numpy as np
    M = np.shape(face_points)[0]
    gaps = find_gaps(face_points)
    points_mod = np.zeros([M+len(gaps),2])
    for k in range(M):
        current = [k,0,[0,1]]
        left = [k-1,0,[0,1]]
        right = [k,1,[0,1]]
        if k not in gaps:
            length_left = length(face_points[left],face_points[current])
            length_right = length(face_points[current],face_points[right])
            len_l = length_left / (length_left+length_right)
            len_r = length_right / (length_left+length_right)
            sens_l = sens_map_wall[k-1,:2] * len_l
            sens_r = sens_map_wall[k,:2] * len_r

            points_mod[k,:2] = sens_l + sens_r

        else:
            pass
    return points_mod
