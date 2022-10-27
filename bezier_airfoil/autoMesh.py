import numpy as np
import matplotlib.pyplot as plt
from bezier_foil import *
import sys

# This function provides an optional command line input to clean up the folder in case of old test cases needing deletion.

retain = ['0', 'constant', 'system', 'data_process.ipynb', 'autoMesh.py', 'output','bezier_foil.py', '.ipynb_checkpoints']

# Define Chord Length for all Other Scaling:
chord_length = 0.3
def store(retain):
    import os
    import shutil
    ignore = ['data_process.ipynb','autoMesh.py','bezier_foil.py','.ipynb_checkpoints']
    copy = ['0','system','constant','dynamicCode']
    dirs = os.listdir()
    delete = []
    for item in dirs:
        if item[0:3] == 'pro':
            delete.append(item)
    target = input("input case name: ")
    casefile = "/home/james/Documents/research/completed_cases/coanda_airfoils/{}/".format(target)
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
    print(retain)
    import shutil
    import os
    dirs = os.listdir()
    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)
        else:
            continue

#---------- Initialization of Coanda Definition --------------#

args = sys.argv

def arg_handle(args):
    # Default Values
    scale = 1
    Rc = 0.14 * .0254       # Coanda cylinder radius
    te = 0.007 * .0254      # R/te = 20
    tu = 0.009 * .0254      # Upper surface thickness at exit

    for i in range(np.size(args)):
        current = args[i].lower()
        if current == '-scale':
            scale = float(args[i+1])
            Rc = Rc * scale
            te = te * scale
            tu = tu * scale
        elif current == '-save':
            store(retain)
        elif current == '-clean':
            cleanup(retain)
    return Rc, te, tu

r,h,t = arg_handle(args)

#---------- Control Variables handled by this block ----------#

control_points_init = np.array([[0,0],[0,.05],[.25,.05],[.5,.07],[.75,.03],[1,0]])
file = '/home/james/Documents/research/cfd/airfoils/naca0015-il.csv'
bezfoil, cpU, iters = foil_opt(control_points_init, file, chord_length, 5e-7,step = 2,debug=True)
# Command line argument handler:
cpL = cpU * [1,-1]
m = 101

upper, lower, coanda = coanda_foil(r,t,h,cpU)

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


for i in range(1,np.shape(upper)[0]-1):
    # pback = np.vstack([pback,cpU[i,:]])
    cpu_string_b += np.array2string(np.hstack([upper[i,:],bk])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'
    cpu_string_f += np.array2string(np.hstack([upper[i,:],fr])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'

for j in range(1,np.shape(lower)[0]-1):
    # pback = np.vstack([pback,cpL[j,:]])
    cpl_string_b += np.array2string(np.hstack([lower[j,:],bk])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'
    cpl_string_f += np.array2string(np.hstack([lower[j,:],fr])).replace('[','\t\t\t\t\t\t\t(').replace(']',')') + '\n'




# Finally: Define the blocking and grading parameters!

blocks_x_L = 100
blocks_y_L = 80
blocks_x_R = 80
blocks_y_R = 80
blocks_y_co = 5
blocks_x_flat = 25
grade_x_L = 10
grade_y = 1500

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
    '\t\thex ( 4 28  8  0  5 29  9  1) ({} {} 1) edgeGrading (50 1 1 50 {} {} {} {} 1 1 1 1) // 0\n'.format(blocks_x_L, blocks_y_L, grade_y, grade_y, grade_y, grade_y),
    '\t\thex (28 26 10  8 29 27 11  9) ({} {} 1) simpleGrading ({}  {}  1) // 1\n'.format(blocks_x_flat, blocks_y_R, 1/grade_x_L, grade_y),
    '\t\thex (26 34 14 10 27 35 15 11) ({} {} 1) simpleGrading ( 1  {}  1) // 2\n'.format(blocks_x_R, blocks_y_R, grade_y),
    '\t\thex (34 36  2 14 35 37  3 15) ({} {} 1) simpleGrading ({}  {}  1) // 3\n'.format(blocks_x_flat, blocks_y_L, grade_x_L, grade_y),
    '\t\thex (36  4  0  2 37  5  1  3) ({} {} 1) edgeGrading (.02 1 1 .02 {} {} {} {} 1 1 1 1) // 4\n'.format(blocks_x_L, blocks_y_L, grade_y, grade_y, grade_y, grade_y),
    '\t\thex (24 32 34 26 25 33 35 27) ({} {} 1) simpleGrading ( 1  {}  1) // 5\n'.format(blocks_x_R, blocks_y_co, 1),
    '\t\thex (22 30 32 24 23 31 33 25) ({} {} 1) simpleGrading ( 1  {}  1) // 6\n'.format(blocks_x_R, blocks_y_co, 1),
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
