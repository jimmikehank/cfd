import numpy as np
import matplotlib.pyplot as plt
from bezier_foil import *
import sys
from os import path

# This function provides an optional command line input to clean up the folder in case of old test cases needing deletion.

retain = ['0', 'constant', 'system', 'data_process.ipynb', 'autoMesh.py', 'output','bezier_foil.py', '.ipynb_checkpoints','autoCFD.py','airfoil_comparison.png']
U_filename = './0/U'
T_filename = './0/T'

# Define Chord Length for all Other Scaling:
chord_length = 1
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
    casefile = "/media/james/Data/james/completed_cases/bezier_airfoils/temps/{}/".format(target)
    print(casefile)
    if os.path.exists(casefile):
        existing = os.listdir(casefile)
        print(existing)
    else:
        os.mkdir(casefile)
        print('make')
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

def change_line_temp(T_filename,temp):
    with open(T_filename,'r') as f:
        test_lines = f.readlines()
        f.close()
    new_line = '\t\t\t\tvalue\t\t\t\t\t\tuniform {};\t\t\t\t\t\t // temperature of wing\n'.format(temp)
    test_lines[59] = new_line
    with open(T_filename,'w') as g:
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
        elif arg[i] == '-temp':
            change_line_temp(T_filename,arg[i+1])
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
symmetry = True
file_exists = path.isfile('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
if file_exists:
    cpU = np.loadtxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel))
    cpL = np.loadtxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel))
    cpU = cpU * chord_length/(cpU[-1,0]-cpU[0,0])
    cpL = cpL * chord_length/(cpL[-1,0]-cpL[0,0])
    bezfoil = init_bezfoil(m,cpU,control_points_lower=cpL,symmetric=False)
else:
    file = '/home/james/Documents/research/cfd/airfoils/{}-il.csv'.format(airfoil_sel)
    bezfoil, cpU, cpL, iters = foil_opt(control_points_init, file, chord_length, 1e-6,m=m,step = 2,debug=True,control_points_lower=cpl_init,sym=symmetry)
    np.savetxt('{}/control_points/{}_cpu.txt'.format(airfoil_dir,airfoil_sel),cpU)
    np.savetxt('{}/control_points/{}_cpl.txt'.format(airfoil_dir,airfoil_sel),cpL)
    print(cpU[-1,0]-cpU[0,0])

# Command line argument handler:

upper = bezfoil[0:m,:]
upper = upper[::-1,:]
lower = bezfoil[m:,:]
#-------------------------------------------------------------#

# blockMeshDefaults:
scale = 1

# Other important values:
span = 1
farfield = 25 * chord_length
trailing = 6 * chord_length
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
pback[1,:] = np.array([0,0])
pback[2,:] = np.array([trailing,bU])
pback[3,:] = np.array([-radius,0])
pback[4,:] = np.array([trailing,bD])
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

blocks_x = 40
blocks_y = 40
grade_x = 1500
egrade_x = 20
egrade_o = 1
grade_y = 200

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
    '\t\thex ( 0  2  4  6  1  3  5  7 ) ({} {} 1) edgeGrading (((.5 .5 {}) (.5 .5 {})) ((.5 .5 {}) (.5 .5 {}))   ((.5 .5 {}) (.5 .5 {})) ((.5 .5 {}) (.5 .5 {})) {} {} {} {} 1 1 1 1) // 0\n'.format(blocks_x*2, blocks_y, egrade_x, 1/egrade_x, egrade_o, 1/egrade_o, egrade_o, 1/egrade_o, egrade_x, 1/egrade_x, 2*grade_y, grade_y, grade_y, 2*grade_y),
    '\t\thex ( 6  8 16  0  7  9 17  1 ) ({} {} 1) edgeGrading (((.5 .5 {}) (.5 .5 {})) ((.5 .5 {}) (.5 .5 {}))   ((.5 .5 {}) (.5 .5 {})) ((.5 .5 {}) (.5 .5 {})) {} {} {} {} 1 1 1 1) // 1\n'.format(blocks_x*2, blocks_y, egrade_o, 1/egrade_o, egrade_x, 1/egrade_x, egrade_x, 1/egrade_x, egrade_o, 1/egrade_o, 0.5/grade_y, 1/grade_y, 1/grade_y, 0.5/grade_y),
    '\t\thex ( 8 10 12 16  9 11 13 17 ) ({} {} 1) edgeGrading ( 1 {} {} 1 {} {} {} {} 1 1 1 1 ) // 2\n'.format(blocks_x, blocks_y, grade_x, grade_x, 1/grade_y, 1/grade_y, 1/grade_y, 1/grade_y),
    '\t\thex ( 2 12 14  4  3 13 15  5 ) ({} {} 1) edgeGrading ( {} 1 1 {} {} {} {} {} 1 1 1 1 ) // 3\n'.format(blocks_x, blocks_y, grade_x, grade_x, grade_y, grade_y, grade_y, grade_y),
    ');\n\n'
]

edges = [
    'edges\n(\n',
    '\t\tarc  4  6  {}  (0 0 1)\n'.format(theta),
    '\t\tarc  5  7  {}  (0 0 1)\n'.format(theta),
    '\t\tarc  6  8  {}  (0 0 1)\n'.format(theta),
    '\t\tarc  7  9  {}  (0 0 1)\n'.format(theta),
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
