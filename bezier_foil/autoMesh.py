import numpy as np
import matplotlib.pyplot as plt
from bezier_foil import *
import sys
from os import path

# This function provides an optional command line input to clean up the folder in case of old test cases needing deletion.

retain = ['0', 'constant', 'system', 'data_process.ipynb', 'autoMesh.py', 'output','bezier_foil.py', '.ipynb_checkpoints','autoCFD.py','airfoil_comparison.png']
U_filename = './0/U'

# Define Chord Length for all Other Scaling:
chord_length = 0.3
def store(retain, target):
    import os
    import shutil
    ignore = ['data_process.ipynb','autoMesh.py','bezier_foil.py','.ipynb_checkpoints']
    copy = ['0','system','constant','dynamicCode']
    dirs = os.listdir()
    delete = []
    for item in dirs:
        if item[0:3] == 'pro':
            delete.append(item)
    casefile = "/home/james/Documents/research/completed_cases/bezier_airfoils/{}/".format(target)
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

def change_line(U_filename,aoa):
    with open(U_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line = 'alpha\t\t\t\t\t\t{};\t\t\t\t\t\t // AoA in Degrees\n'.format(aoa)
    test_lines[20] = new_line
    with open(U_filename,'w') as g:
        g.writelines(test_lines)
        g.close()

#---------- Initialization of Coanda Definition --------------#

args = sys.argv

def arg_handle(arg):
    # Default Values
    for i in range(len(arg)):
        if arg[i] == "-debug":
            print("Debug mode not implemented")
        elif arg[i] == '-clean':
            cleanup(retain)
        elif arg[i] =='-aoa':
            change_line(U_filename, arg[i+1])
        elif arg[i] == '-save':
            try:
                store(retain,arg[i+1])
            except:
                "No target filename provided"

    return 0
def check_args(arg):
    for i in range(len(arg)):
        if arg[i] == '-aoa':
            aoa = float(arg[i+1])
            return aoa
        else:
            pass
    return 0

argout = arg_handle(args)
aoa = check_args(args)

#---------- Control Variables handled by this block ----------#
m = 101
control_points_init = np.array([[0,0],[0,.05],[.25,.05],[.35,.06],[.5,.07],[.75,.03],[1,0]])
# control_points_init = np.array([[0,0],[0,.05],[.75,.03],[1,0]])
cpl_init = control_points_init * np.array([1,-1])
airfoil_dir = '/home/james/Documents/research/cfd/airfoils/'
airfoil_sel = 'rae2822'
symmetry = False
file_exists = path.isfile('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
if file_exists:
    cpu = np.loadtxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
    cpl = np.loadtxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel))
    bezfoil = init_bezfoil(m,cpu,control_points_lower=cpl,symmetric=False)
else:
    file = '/home/james/Documents/research/cfd/airfoils/{}-il.csv'.format(airfoil_sel)
    bezfoil, cpU, cpL, iters = foil_opt(control_points_init, file, chord_length, 1e-6,m=m,step = 2,debug=True,control_points_lower=cpl_init,sym=symmetry)
    np.savetxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel),cpU)
    np.savetxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel),cpL)

# Command line argument handler:

upper = bezfoil[0:m,:]
upper = upper[::-1,:]
lower = bezfoil[m:,:]
#-------------------------------------------------------------#

# blockMeshDefaults:
scale = 1

# Other important values:
span = 1

# Bounding Box:

bU =  10 * chord_length # Upper bound
bD = -10 * chord_length # Lower bound
bL = -10 * chord_length # Left bound
bR =  10 * chord_length # Right bound
fr =  span/2            # Front bound
bk = -span/2            # Back bound


# print(fil,p0)
pback = np.zeros([9,2])
eps = 1e-6
back_deflection = bR*np.sin(np.pi*aoa/180)

# Define Point Matrix:
pback[0,:] = np.array([-chord_length,0])
pback[1,:] = np.array([0,0])
pback[2,:] = np.array([0,bU])
pback[3,:] = np.array([bL,0])
pback[4,:] = np.array([0,bD])
pback[5,:] = np.array([bR,bL])
pback[6,:] = np.array([bR,back_deflection])
pback[7,:] = np.array([bR,bU])
pback[8,:] = np.array([0,0])

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

blocks_x = 60
blocks_y = 40
grade_x = 200
egrade_x = 5
grade_y = 1000

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

blocks = [
    'blocks\n(\n',
    '\t\thex ( 0  2  4  6  1  3  5  7 ) ({} {} 1) edgeGrading ({}   1   1  {} {} {} {} {} 1 1 1 1) // 0\n'.format(blocks_x*2, blocks_y, egrade_x, egrade_x, grade_y, grade_y, grade_y, grade_y),
    '\t\thex ( 6  8 16  0  7  9 17  1 ) ({} {} 1) edgeGrading ( 1  {}  {}   1 {} {} {} {} 1 1 1 1) // 1\n'.format(blocks_x*2, blocks_y, egrade_x, egrade_x, 1/grade_y, 1/grade_y, 1/grade_y, 1/grade_y),
    '\t\thex ( 8 10 12 16  9 11 13 17 ) ({} {} 1) simpleGrading ( {} {} 1) // 1\n'.format(blocks_x, blocks_y, grade_x, 1/grade_y),
    '\t\thex ( 2 12 14  4  3 13 15  5 ) ({} {} 1) simpleGrading ( {} {} 1) // 1\n'.format(blocks_x, blocks_y, grade_x, grade_y),
    ');\n\n'
]

edges = [
    'edges\n(\n',
    '\t\tarc  4  6  90.0  (0 0 1)\n',
    '\t\tarc  5  7  90.0  (0 0 1)\n',
    '\t\tarc  6  8  90.0  (0 0 1)\n',
    '\t\tarc  7  9  90.0  (0 0 1)\n',
    '\t\tspline 0  2 (\n{}\n\t\t\t\t\t\t)\n'.format(cpu_string_b),
    '\t\tspline 1  3 (\n{}\n\t\t\t\t\t\t)\n'.format(cpu_string_f),
    '\t\tspline 0 16 (\n{}\n\t\t\t\t\t\t)\n'.format(cpl_string_b),
    '\t\tspline 1 17 (\n{}\n\t\t\t\t\t\t)\n'.format(cpl_string_f),
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
