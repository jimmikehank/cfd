import numpy as np
import matplotlib.pyplot as plt
from bezier_foil import *
import sys
from os import path
import argparse
import time

# This function provides an optional command line input to clean up the folder in case of old test cases needing deletion.

retain = ['0', 'constant', 'system', 'data_process.ipynb', 'autoMesh.py', 'output','bezier_foil.py', '.ipynb_checkpoints','autoCFD.py','airfoil_comparison.png','processing.py']
U_filename = './0/U'
# Define Chord Length for all Other Scaling:
chord_length = 0.3

parser = argparse.ArgumentParser()
parser.add_argument('--iter', default = 0, type = int, help = 'Iteration for selection of optimization control points')
parser.add_argument('--airfoil', default = 'naca0015', type = str, help = 'Selection of initial airfoil')
parser.add_argument('--runName', default = 'fail_check', type = str, help = 'Name of target directory for control points')
parser.add_argument('--aoa', default = 0., type = float, help = 'Specified angle of attack in degrees, \nDefault: 0.0')
parser.add_argument('--clean', default = False, type = bool, help = 'Run clean function. \nDefault: False')
parser.add_argument('--store', default = False, type = bool, help = 'Run store function. \nDefault: False')
parser.add_argument('--mdot', default = 0.0, type = float, help = 'Mass flow for coanda in kg/s')
parser.add_argument('--meanFlow', default = 25.0, type = float, help = 'Mean flow speed for simulation in m/s \n Default: 25.0')
parser.add_argument('--control', default = 0, type = int, help = 'Selected variable to step in for forward finite differencing. \nValues:\n0: None\n1:Coanda Radius\n2:Coanda Opening Height\n3: Knife Edge Thickness')
parser.add_argument('--delta', default = 1e-5, type = float, help = 'Size of step for forward differencing\nDefault: 1e-5')
parser.add_argument('--pressUp', default = 100001, type = float, help = 'Pressure in upper plenum')
parser.add_argument('--pressLo', default = 100001, type = float, help = 'Pressure in lower plenum')
args = parser.parse_args()

clean = args.clean
airfoil = args.airfoil.lower()
aoa = args.aoa
meanflow = args.meanFlow
savebool = args.store
controls = args.control
iters = args.iter
delta = args.delta
pupper = args.pressUp
plower = args.pressLo
runName = args.runName

def store(retain,name):
    import os
    import shutil
    import time
    # print(retain)
    print("Saving run in {}".format(name))
    time.sleep(2)
    ignore = ['data_process.ipynb','autoMesh.py','bezier_foil.py','.ipynb_checkpoints','processing.py']
    copy = ['0','system','constant','dynamicCode']
    dirs = os.listdir()
    delete = []
    for item in dirs:
        if item[0:9] == 'processor':
            delete.append(item)
    target = name
    casefile = "/home/james/Documents/research/completed_cases/pl_redo/{}/".format(target)
    if os.path.exists(casefile):
        existing = os.listdir(casefile)
    else:
        os.system("mkdir -p {}".format(casefile))
        existing = []
    for item in dirs:
        if item in delete:
            shutil.rmtree(item)
        elif item not in retain and item not in ignore and item not in existing:
            shutil.move(item,casefile)
        elif item in copy:
            shutil.copytree(item,casefile+item)


def cleanup(retain):
    import shutil
    import os
    dirs = os.listdir()
    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)
        else:
            continue

def change_line_aoa(U_filename,aoa):
    with open(U_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line = 'alpha\t\t\t\t\t\t{};\t\t\t\t // Angle of Attack\n'.format(aoa)
    test_lines[21] = new_line
    with open(U_filename,'w') as g:
        g.writelines(test_lines)
        g.close()

def change_line_meanflow(U_filename,meanflow):
    with open(U_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line = 'U\t\t\t\t\t\t\t\t{};\t\t\t\t\t\t // Mean Flow Speed\n'.format(meanflow)
    test_lines[22] = new_line
    with open(U_filename,'w') as g:
        g.writelines(test_lines)
        g.close()

def change_line_pressure(P_filename,pupper,plower):
    with open(P_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line_u = 'upper_pressure uniform {};\t\t\t\t\t\t // Upper Plenum Pressure\n'.format(pupper)
    new_line_l = 'lower_pressure uniform {};\t\t\t\t\t\t // Lower Plenum Pressure\n'.format(plower)
    test_lines[22] = new_line_u
    test_lines[23] = new_line_l
    with open(P_filename,'w') as g:
        g.writelines(test_lines)
        g.close()
#---------- Initialization of Coanda Definition --------------#


Ufile = './0/U'
Pfile = './0/p'
if clean == True:
    cleanup(retain)
if savebool == True:
    store(retain,runName)

change_line_aoa(Ufile,aoa)
change_line_meanflow(Ufile,meanflow)
change_line_pressure(Pfile, pupper, plower)
print("Changed pressure to {}".format(pupper))



def loadnumbers(iters):
    filename = '/home/james/Documents/research/completed_cases/coanda_opt/{}.txt'.format(iters)
    data = np.loadtxt(filename)
    Rc = data[0]
    te = data[1]
    tu = data[2]
    return Rc, te, tu

def coanda_set(iters,control,delta):
    # Default Values
    if iters == 0:
        Rc = 0.12 * .0254       # Coanda cylinder radius
        te = 0.009 * .0254      # R/te = 20
        tu = 0.012 * .0254      # Upper surface thickness at exit
    else:
        Rc, te, tu = loadnumbers(iters)

    if control == 0:
        return Rc, te, tu
    elif control == 1:
        Rc = Rc + delta
    elif control == 2:
        te = te + delta
    elif control == 3:
        tu = tu + delta
    return Rc, te, tu

r,h,t = coanda_set(iters,controls,delta)

print("Coanda Dimensions: R: {}, h: {}, t: {}".format(r,h,t))
time.sleep(1)
#---------- Control Variables handled by this block ----------#
#---------- Control Variables handled by this block ----------#
m = 101
control_points_init = np.array([[0,0],[0,.05],[.25,.05],[.35,.06],[.5,.07],[.75,.03],[1,0]])
# control_points_init = np.array([[0,0],[0,.05],[.75,.03],[1,0]])
cpl_init = control_points_init * np.array([1,-1])
airfoil_dir = '/home/james/Documents/research/cfd/airfoils/'
airfoil_sel = airfoil
symmetry = False
file_exists = path.isfile('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
if file_exists:
    cpu = np.loadtxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
    cpl = np.loadtxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel))
    bezfoil = init_bezfoil(m,cpu,control_points_lower=cpl,symmetric=False)
else:
    file = '/home/james/Documents/research/cfd/airfoils/{}-il.csv'.format(airfoil_sel)
    bezfoil, cpu, cpl, iters = foil_opt(control_points_init, file, chord_length, 1e-6,m=m,step = 2,debug=True,control_points_lower=cpl_init,sym=symmetry)
    np.savetxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel),cpu)
    np.savetxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel),cpl)

# Command line argument handler:

upper, lower, coanda = coanda_foil(r,t,h,cpu,cpl)

#-------------------------------------------------------------#

# blockMeshDefaults:
scale = 1

# Other important values:
slot_width = .156

# Bounding Box:

bU =  10 * chord_length # Upper bound
bD = -10 * chord_length # Lower bound
bL = -10 * chord_length # Left bound
bR =  10 * chord_length # Right bound
fr =  slot_width/2                 # Front bound
bk = -slot_width/2                 # Back bound

# print(fil,p0)
pback = np.zeros([10,2])
eps = 1e-6

# New coordinates for Plenum:
Rc = r
te = h
tu = t
le = -r
pln = 0.04              # Plenum length
pcn = (Rc - te)*1.5
hcn = (Rc-te)/2
beta = np.arctan(pcn/hcn) * 180 / np.pi
ai = 180 - 2*beta
# print(fil,p0)

# Define Point Matrix:
pback[0,:] = np.array([bL, 0])
pback[1,:] = np.array([lower[-1,0], bD])
pback[2,:] = upper[0,:]
pback[3,:] = upper[-1,:]
pback[4,:] = np.array([upper[-1,0],bU])
pback[5,:] = np.array([-r+3*eps,bU])
pback[6,:] = np.array([bR-r,0])
pback[7,:] = np.array([-r+3*eps,bD])
pback[8,:] = np.array([0,0])
pback[9,:] = lower[-1,:]

for i in range(np.shape(coanda)[0]):
    pback = np.vstack([pback,coanda[i,:]])
cpu_string_b = ""
cpu_string_f = ""
cpl_string_b = ""
cpl_string_f = ""

appstack = np.zeros([12,2])
pback = np.vstack([pback,appstack])

pback[19,:] = np.array([-pln+le, Rc+te])
pback[20,:] = np.array([-pln+le, te])
pback[21,:] = np.array([-2*pcn+le,te])
pback[22,:] = np.array([-2*pcn+le,Rc+te])
pback[23,:] = np.array([-pcn+le,hcn+te])
pback[24,:] = np.array([-pcn+le,Rc+te])

pback[25,:] = np.array([-pln+le, -(Rc+te)])
pback[26,:] = np.array([-pln+le, -te])
pback[27,:] = np.array([-2*pcn+le,-te])
pback[28,:] = np.array([-2*pcn+le,-(Rc+te)])
pback[29,:] = np.array([-pcn+le,-(hcn+te)])
pback[30,:] = np.array([-pcn+le,-(Rc+te)])

for i in range(1,np.shape(upper)[0]-1):
    # pback = np.vstack([pback,cpU[i,:]])
    cpu_string_b += np.array2string(np.hstack([upper[i,:],bk])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'
    cpu_string_f += np.array2string(np.hstack([upper[i,:],fr])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'

for j in range(1,np.shape(lower)[0]-1):
    # pback = np.vstack([pback,cpL[j,:]])
    cpl_string_b += np.array2string(np.hstack([lower[j,:],bk])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'
    cpl_string_f += np.array2string(np.hstack([lower[j,:],fr])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'




# Finally: Define the blocking and grading parameters!

blocks_x_L = 90
blocks_y_L = 110
blocks_x_R = 150
blocks_y_R = blocks_y_L
blocks_y_co = 5
blocks_x_flat = 20
blocks_x_in = 10
blocks_y_in = blocks_y_co
grade_x_L = 2
egrade_x = 10
egrade_o = 10
grade_y = 1600

header = [
    '/*---------------------------------*- C++ -*-----------------------------------*/\n',
    '//    //////////////      ///////////////\n',
    '//    //////////////      //////////////\n',
    '//         ////           ////             /////       ////  //  /// //// ////\n',
    '//         ////     ////  ///////////   ///    ///  ///    ///   ////  ///  ///\n',
    '//  ////   ////           ///////////   ///    ///  ///    ///   ///   ///  ///\n',
    '//   //// ////            ////          ///    ///  ///    ///   ///   ///  ///\n',
    '//     /////              ////            /////       ///// ///  ///   ///  ///\n',
    '/*---------------------------------*- C++ -*-----------------------------------*/\n\n'
]

opener = [
    'FoamFile',
    '\n',
    '{',
    '\n',
    '\t\tversion\t\t\t2.0;',
    '\n',
    '\t\tformat\t\t\tascii;',
    '\n',
    '\t\tclass\t\t\t\tdictionary;',
    '\n',
    '\t\tobject\t\t\tblockMeshDict;',
    '\n',
    '}\n\n',
    '// ***************************************************************************** //\n\n'
]

scale = 'scale\t\t\t\t1;\n\n'

pointsdef = []

def makepoints(points_back,back,front):
    pointsdef = ['vertices\n','(\n']
    for i in range(np.shape(pback)[0]*2):
        if i % 2 == 0:
            pointsdef.append('\t\t({} {} {})\t//\t{}\n'.format(points_back[int(i/2),0],points_back[int(i/2),1],back,i))

        else:
            pointsdef.append('\t\t({} {} {})\t//\t{}\n'.format(points_back[int((i-1)/2),0],points_back[int((i-1)/2),1],front,i))

    pointsdef.append(');\n\n')
    return pointsdef

# Block order is inside coanda plenum 0, 1, 2 then 3 is the block directly below the surface, the rest follow counter-clockwise


blocks = [
    'blocks\n(\n',
    '\t\thex ( 4 28  8  0  5 29  9  1) ({} {} 1) edgeGrading (((.5 .5 {}) (.5 .5 {})) {}  {} ((.5 .5 {}) (.5 .5 {})) {} {} {} {} 1 1 1 1) // 0\n'.format(blocks_x_L, blocks_y_L, egrade_x, 1/egrade_x, 1/egrade_o, 1/egrade_o, egrade_x, 1/egrade_x, grade_y, grade_y, grade_y, grade_y),
    '\t\thex (28 26 10  8 29 27 11  9) ({} {} 1) simpleGrading ({}  {}  1) // 1\n'.format(blocks_x_flat, blocks_y_R, 1/grade_x_L, grade_y),
    '\t\thex (26 34 14 10 27 35 15 11) ({} {} 1) simpleGrading ( 1  {}  1) // 2\n'.format(blocks_x_R, blocks_y_R, grade_y),
    '\t\thex (34 36  2 14 35 37  3 15) ({} {} 1) simpleGrading ({}  {}  1) // 3\n'.format(blocks_x_flat, blocks_y_L, grade_x_L, grade_y),
    '\t\thex (36  4  0  2 37  5  1  3) ({} {} 1) edgeGrading (((.5 .5 {}) (.5 .5 {})) {}  {} ((.5 .5 {}) (.5 .5 {})) {} {} {} {} 1 1 1 1) // 4\n'.format(blocks_x_L, blocks_y_L, egrade_x, 1/egrade_x, egrade_o, egrade_o, egrade_x, 1/egrade_x, grade_y, grade_y, grade_y, grade_y),
    '\t\thex (24 32 34 26 25 33 35 27) ({} {} 1) simpleGrading ( 1  {}  1) // 5\n'.format(blocks_x_R, blocks_y_co-2, 2),
    '\t\thex (22 30 32 24 23 31 33 25) ({} {} 1) simpleGrading ( 1  {}  1) // 6\n'.format(blocks_x_R, blocks_y_co, 1),
    '\t\thex (40 42 44 38 41 43 45 39) ({} {} 1) simplegrading ( 1   1  1) // 7\n'.format(blocks_x_in+20,blocks_y_in),
    '\t\thex (42 46 48 44 43 47 49 45) ({} {} 1) simpleGrading ( 1   1  1) // 8\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (46 22 24 48 47 23 25 49) ({} {} 1) simpleGrading ( 1   1  1) // 9\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (50 56 54 52 51 57 55 53) ({} {} 1) simpleGrading ( 1   1  1) //10\n'.format(blocks_x_in+20,blocks_y_in),
    '\t\thex (56 60 58 54 57 61 59 55) ({} {} 1) simpleGrading ( 1   1  1) //11\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (60 32 30 58 61 33 31 59) ({} {} 1) simpleGrading ( 1   1  1) //12\n'.format(blocks_x_in,blocks_y_in),
    ');\n\n'
]

edges = [
    'edges\n(\n',
    '\t\tarc  8  0  90.0  (0 0 1)\n',
    '\t\tarc  9  1  90.0  (0 0 1)\n',
    '\t\tarc  0  2  90.0  (0 0 1)\n',
    '\t\tarc  1  3  90.0  (0 0 1)\n',
    '\t\tarc 14 10 180.0  (0 0 1)\n',
    '\t\tarc 15 11 180.0  (0 0 1)\n',
    '\t\tarc 30 22 180.0  (0 0 1)\n',
    '\t\tarc 31 23 180.0  (0 0 1)\n',
    '\t\tarc 32 24 180.0  (0 0 1)\n',
    '\t\tarc 33 25 180.0  (0 0 1)\n',
    '\t\tarc 34 26 180.0  (0 0 1)\n',
    '\t\tarc 35 27 180.0  (0 0 1)\n',
    '\t\tspline 4 28 (\n{}\n\t\t\t\t\t\t)\n'.format(cpu_string_b),
    '\t\tspline 4 36 (\n{}\n\t\t\t\t\t\t)\n'.format(cpl_string_b),
    '\t\tspline 5 29 (\n{}\n\t\t\t\t\t\t)\n'.format(cpu_string_f),
    '\t\tspline 5 37 (\n{}\n\t\t\t\t\t\t)\n'.format(cpl_string_f),
    '\t\tarc  42 46    {} (0 0  1)\n'.format(ai),
    '\t\tarc  46 22    {} (0 0 -1)\n'.format(ai),
    '\t\tarc  43 47    {} (0 0  1)\n'.format(ai),
    '\t\tarc  47 23    {} (0 0 -1)\n'.format(ai),

    '\t\tarc  54 58    {} (0 0 -1)\n'.format(ai),
    '\t\tarc  58 30    {} (0 0  1)\n'.format(ai),
    '\t\tarc  55 59    {} (0 0 -1)\n'.format(ai),
    '\t\tarc  59 31    {} (0 0  1)\n'.format(ai),
    ');\n\n'
]

pointsdef = makepoints(pback,bk,fr)

with open('./system/blockMeshDict','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(scale)
    f.writelines(pointsdef)
    f.writelines(blocks)
    f.writelines(edges)
    with open('./system/faces','r') as g:
        f.writelines(g.readlines())
