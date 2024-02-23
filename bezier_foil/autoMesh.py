import numpy as np
import matplotlib.pyplot as plt
from bezier_foil import *
import sys
from os import path
import argparse

# This function provides an optional command line input to clean up the folder in case of old test cases needing deletion.

retain = ['0', 'constant', 'system', 'data_process.ipynb', 'autoMesh.py', 'output','bezier_foil.py', 'processing.py', '.ipynb_checkpoints','autoCFD.py','airfoil_comparison.png','james_test.py','FFD']
file_delete = []
U_filename = './0/U'
T_filename = './0/T'

parser = argparse.ArgumentParser()
parser.add_argument('--iter', default = 0, type = int, help = 'Iteration for selection of optimization control points')
parser.add_argument('--airfoil', default = 'NACA0012', type = str, help = 'Selection of initial airfoil')
parser.add_argument('--runName', default = 'fail_check', type = str, help = 'Name of target directory for control points')
parser.add_argument('--aoa', default = 0., type = float, help = 'Specified angle of attack in degrees, \nDefault: 0.0')
parser.add_argument('--clean', default = False, type = bool, help = 'Run clean function. \nDefault: False')
parser.add_argument('--store', default = False, type = bool, help = 'Run store function. \nDefault: False')
parser.add_argument('--meanFlow', default = 25.0, type = float, help = 'Mean flow speed for simulation\nDefault: 25 m/s')
parser.add_argument('--finiteTE', default = False, type = bool, help = 'Enables finite TE for meshing\nDefault: False')

args = parser.parse_args()
iteration = args.iter
airfoil_sel = args.airfoil.lower()
runName = args.runName
aoa = args.aoa
clean_bool = args.clean
store_bool = args.store
meanflow = args.meanFlow
finiteTE = args.finiteTE

# Define Chord Length for all Other Scaling:
chord_length = 1.0
# TE Thickness / Chord Length
te_c_ratio = 0.002
# Spatial Resolution Multiplier
spatial_mul = 1

def store(retain, target):
    import os
    import shutil
    ignore = ['data_process.ipynb','autoMesh.py','bezier_foil.py','.ipynb_checkpoints']
    copy = ['0','system','constant','dynamicCode']
    dirs = os.listdir()
    delete = []
    for item in dirs:
        if item[0:9] == 'processor':
            delete.append(item)
    casefile = target
    if os.path.exists(casefile):
        existing = os.listdir(casefile)
    else:
        os.mkdir(casefile)
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
        if item not in retain and os.path.isdir(item):
            shutil.rmtree(item)
        elif item not in retain and os.path.isfile(item):
            os.remove(item)
        else:
            continue

def change_line(U_filename,aoa):
    with open(U_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line = 'alpha\t\t\t\t\t\t{};\t\t\t\t\t\t // AoA in Degrees\n'.format(aoa)
    test_lines[20] = new_line
    with open(U_filename,'w') as g:
        g.writelines(test_lines)
        g.close()

def change_line_meanflow(U_filename,meanflow):
    with open(U_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line = 'U\t\t\t\t\t\t{};\t\t\t\t\t\t // AoA in Degrees\n'.format(meanflow)
    test_lines[21] = new_line
    with open(U_filename,'w') as g:
        g.writelines(test_lines)
        g.close()

#---------- Initialization of Argument Variables -------------#

m = 101
airfoil_dir = '/home/james/Documents/research/cfd/airfoils/'
change_line(U_filename,aoa)
change_line_meanflow(U_filename,meanflow)

if clean_bool:
    cleanup(retain)
if store_bool:
    target = '/{}/'.format(runName)
    store(retain,target)

# First iteration of optimization loop, run initial shape match to airfoil for bezier curves.
if iteration == 0:
    control_points_init = np.array([[0,0],[0,.05],[.25,.05],[.35,.06],[.5,.06],[.75,.03],[1,0]])
    cpl_init = control_points_init * np.array([1,-1])
    symmetry = False
    file_exists = path.isfile('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
    if file_exists:
        cpU = np.loadtxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
        cpL = np.loadtxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel))
        cpU = cpU * chord_length/(cpU[-1,0]-cpU[0,0])
        cpL = cpL * chord_length/(cpL[-1,0]-cpL[0,0])
        bezfoil = init_bezfoil(m,cpU,control_points_lower=cpL,symmetric=symmetry)
    else:
        file = '/home/james/Documents/research/cfd/airfoils/{}-il.csv'.format(airfoil_sel)
        bezfoil, cpU, cpL, iters = foil_opt(control_points_init, file, chord_length, 1e-6, m = m, step = 2, debug = True, control_points_lower = cpl_init, sym = symmetry)
        np.savetxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel),cpU)
        np.savetxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel),cpL)
        print(cpU[-1,0]-cpU[0,0])
    if finiteTE:
        cpU[-1,1] += te_c_ratio * chord_length
        cpL[-1,1] -= te_c_ratio * chord_length
        upper = bezier_curve(cpU, m)
        lower = bezier_curve(cpL, m)
    else:
        upper = bezfoil[0:m,:]
        upper = upper[::-1,:]
        lower = bezfoil[m:,:]

else:
    cpU = np.loadtxt('{}/optimization/{}/{}_cpu.txt'.format(airfoil_dir, runName, iteration-1))
    cpL = np.loadtxt('{}/optimization/{}/{}_cpl.txt'.format(airfoil_dir, runName, iteration-1))

    cpU = cpU * chord_length/(cpU[-1,0]-cpU[0,0])
    cpL = cpL * chord_length/(cpL[-1,0]-cpL[0,0])

    bezfoil = init_bezfoil(m, cpU, control_points_lower = cpL, symmetric = False)
    upper = bezfoil[0:m,:]
    upper = upper[::-1,:]
    lower = bezfoil[m:,:]
#-------------------------------------------------------------#

# blockMeshDefaults:
scale = 1

# Other important values:
span = 1
farfield = 25 * chord_length
trailing = 8 * chord_length

# Automated Calculation, no need to touch this one
theta = 90 + np.arctan(trailing/farfield)*(180/np.pi)

# Bounding Box:
bU =  farfield          # Upper bound
bD = -farfield          # Lower bound
bL = -farfield          # Left bound
bR =  farfield          # Right bound
fr =  span/2            # Front bound
bk = -span/2            # Back bound

radius = np.sqrt(bU**2 + trailing**2)


# print(fil,p0)
pback = np.zeros([9,2])
eps = 1e-6
back_deflection = bR*np.sin(np.pi*aoa/180)

# Define Point Matrix:
pback[0,:] = np.array([-chord_length,0])
pback[1,:] = np.array([0,upper[-1,1]])
pback[2,:] = np.array([trailing,bU])
pback[3,:] = np.array([-radius,0])
pback[4,:] = np.array([trailing,bD])
pback[5,:] = np.array([bR,bL])
pback[6,:] = np.array([bR,back_deflection])
pback[7,:] = np.array([bR,bU])
pback[8,:] = np.array([0,lower[-1,1]])

if finiteTE:
    newpts = np.zeros([2,2])
    newpts[0,:] = np.array([bR,upper[-1,1]*10+back_deflection])
    newpts[1,:] = np.array([bR,lower[-1,1]*10+back_deflection])
    pback = np.vstack([pback,newpts])

cpu_string_b = ''
cpu_string_f = ''
cpl_string_b = ''
cpl_string_f = ''

for i in range(1,np.shape(upper)[0]-1):
    # pback = np.vstack([pback,cpU[i,:]])
    cpu_string_b += np.array2string(np.hstack([upper[i,:],bk])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'
    cpu_string_f += np.array2string(np.hstack([upper[i,:],fr])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'

for j in range(1,np.shape(lower)[0]-1):
    # pback = np.vstack([pback,cpL[j,:]])
    cpl_string_b += np.array2string(np.hstack([lower[j,:],bk])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'
    cpl_string_f += np.array2string(np.hstack([lower[j,:],fr])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'




# Finally: Define the blocking and grading parameters!

blocks_x = 50 * spatial_mul
blocks_x_foil = 50 * spatial_mul
blocks_y = 70 * spatial_mul
blocks_TE = 1 * spatial_mul
grade_x = 1200
egrade_x = 10
egrade_o = 10
grade_y = 800
grade_yo = 25
grmul = 10

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



def makepoints(points_back,back,front):
    pointsdef = ['vertices\n','(\n']
    for i in range(np.shape(pback)[0]*2):
        if i % 2 == 0:
            pointsdef.append('\t\t({} {} {})\t//\t{}\n'.format(points_back[int(i/2),0],points_back[int(i/2),1],back,i))

        else:
            pointsdef.append('\t\t({} {} {})\t//\t{}\n'.format(points_back[int((i-1)/2),0],points_back[int((i-1)/2),1],front,i))

    pointsdef.append(');\n\n')
    return pointsdef

pointsdef = makepoints(pback,bk,fr)

# Block order is inside coanda plenum 0, 1, 2 then 3 is the block directly below the surface, the rest follow counter-clockwise

if finiteTE:
    blocks = [
        'blocks\n(\n',
        '\t\thex ( 0  2  4  6  1  3  5  7 ) ({} {} 1) edgeGrading (((.5 .5 {}) (.5 .5 {}))   {}   {}   ((.5 .5 {}) (.5 .5 {})) {} {} {} {} 1 1 1 1) // 0\n'.format(blocks_x_foil, blocks_y, egrade_x, 1/egrade_x, 1/egrade_o, 1/egrade_o, egrade_x, 1/egrade_x, grmul*grade_y, grade_y, grade_y, grmul*grade_y),
        '\t\thex ( 6  8 16  0  7  9 17  1 ) ({} {} 1) edgeGrading ({}   ((.5 .5 {}) (.5 .5 {}))   ((.5 .5 {}) (.5 .5 {}))   {} {} {} {} {} 1 1 1 1) // 1\n'.format(blocks_x_foil, blocks_y, 1/egrade_o, egrade_x, 1/egrade_x, egrade_x, 1/egrade_x, 1/egrade_o, 1/(grmul*grade_y), 1/grade_y, 1/grade_y, 1/(grmul*grade_y)),
        '\t\thex ( 8 10 20 16  9 11 21 17 ) ({} {} 1) edgeGrading ( 1 {} {} 1 {} {} {} {} 1 1 1 1 ) // 2\n'.format(blocks_x, blocks_y, grade_x, grade_x, 1/grade_y, 1/grade_yo, 1/grade_yo, 1/grade_y),
        '\t\thex (16 20 18  2 17 21 19  3 ) ({} {} 1) simpleGrading ( {} 1 1 ) // 3\n'.format(blocks_x, blocks_TE, grade_x),
        '\t\thex ( 2 18 14  4  3 19 15  5 ) ({} {} 1) edgeGrading ( {} 1 1 {} {} {} {} {} 1 1 1 1 ) // 4\n'.format(blocks_x, blocks_y, grade_x, grade_x, grade_y, grade_yo, grade_yo, grade_y),
        ');\n\n'
    ]
else:
    blocks = [
        'blocks\n(\n',
        '\t\thex ( 0  2  4  6  1  3  5  7 ) ({} {} 1) edgeGrading (((.5 .5 {}) (.5 .5 {}))   {}   {}   ((.5 .5 {}) (.5 .5 {})) {} {} {} {} 1 1 1 1) // 0\n'.format(blocks_x_foil, blocks_y, egrade_x, 1/egrade_x, 1/egrade_o, 1/egrade_o, egrade_x, 1/egrade_x, grmul*grade_y, grade_y, grade_y, grmul*grade_y),
        '\t\thex ( 6  8 16  0  7  9 17  1 ) ({} {} 1) edgeGrading ({}   ((.5 .5 {}) (.5 .5 {}))   ((.5 .5 {}) (.5 .5 {}))   {} {} {} {} {} 1 1 1 1) // 1\n'.format(blocks_x_foil, blocks_y, 1/egrade_o, egrade_x, 1/egrade_x, egrade_x, 1/egrade_x, 1/egrade_o, 1/(grmul*grade_y), 1/grade_y, 1/grade_y, 1/(grmul*grade_y)),
        '\t\thex ( 8 10 12 16  9 11 13 17 ) ({} {} 1) edgeGrading ( 1 {} {} 1 {} {} {} {} 1 1 1 1 ) // 2\n'.format(blocks_x, blocks_y, grade_x, grade_x, 1/grade_y, 1/grade_yo, 1/grade_yo, 1/grade_y),
        '\t\thex ( 2 12 14  4  3 13 15  5 ) ({} {} 1) edgeGrading ( {} 1 1 {} {} {} {} {} 1 1 1 1 ) // 3\n'.format(blocks_x, blocks_y, grade_x, grade_x, grade_y, grade_yo, grade_yo, grade_y),
        ');\n\n'
    ]

thetaRad = theta * np.pi/180
mid_x = -radius * np.cos(thetaRad/2)
mid_y = radius * np.sin(thetaRad/2)
edges = [
    'edges\n(\n',
    '\t\tarc  4  6  ({} {} {})\n'.format(mid_x,mid_y,bk),
    '\t\tarc  5  7  ({} {} {})\n'.format(mid_x,mid_y,fr),
    '\t\tarc  6  8  ({} {} {})\n'.format(mid_x,-mid_y,bk),
    '\t\tarc  7  9  ({} {} {})\n'.format(mid_x,-mid_y,fr),
    '\t\tspline 0  2 (\n{}\n\t\t\t\t\t\t)\n'.format(cpu_string_b),
    '\t\tspline 1  3 (\n{}\n\t\t\t\t\t\t)\n'.format(cpu_string_f),
    '\t\tspline 0 16 (\n{}\n\t\t\t\t\t\t)\n'.format(cpl_string_b),
    '\t\tspline 1 17 (\n{}\n\t\t\t\t\t\t)\n'.format(cpl_string_f),
    ');\n\n'
]

pointsdef = makepoints(pback,bk,fr)

if finiteTE:
    facefile = './system/facesTE'
else:
    facefile = './system/faces'

with open('./system/blockMeshDict','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(scale)
    f.writelines(pointsdef)
    f.writelines(blocks)
    f.writelines(edges)
    with open(facefile, 'r') as g:
        f.writelines(g.readlines())
