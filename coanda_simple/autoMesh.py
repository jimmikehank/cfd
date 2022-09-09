import numpy as np
import matplotlib.pyplot as plt
import sys



#---------- Control Variables handled by this block ----------#

retain = ['0', 'constant', 'system', 'data_process.ipynb', 'autoMesh.py', 'output', '.ipynb_checkpoints']

# Command line argument handler:

args = sys.argv

# This function provides an optional command line input to clean up the folder in case of old test cases needing deletion.


def cleanup(retain):
    import shutil
    import os
    dirs = os.listdir()
    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)
        else:
            continue

# This function parses input arguments to generate mesh allowing external programs to iterate through different mesh configurations

def arg_handle(args):
    # Default Values
    Rc = 0.003   # Coanda cylinder radius
    te = 0.00015  # R/te = 20
    tu = 0.00018  # Upper surface thickness at exit

    for i in range(np.size(args)):
        current = args[i].lower()
        if current == '-rc':
            Rc = float(args[i+1])
        elif current == '-te':
            te = float(args[i+1])
        elif current == '-tu':
            tu = float(args[i+1])
        elif current == '-clean':
            cleanup(retain)
    return Rc, te, tu

Rc, te, tu = arg_handle(args)



#-------------------------------------------------------------#

print('\nCoanda Surface Radius = {}\nSlot Exit Width = {}\nExit Surface Thickness = {}\n'.format(Rc, te, tu))

# blockMeshDefaults
scale = 1


# Bounding Box:

bU =  1.0 # Upper bound
bD = -1.0 # Lower bound
bL = -1.0 # Left bound
bR =  1.0 # Right bound
fr =  0.5 # Front bound
bk = -0.5 # Back bound

# Mesh expansion for coanda cylinder

rL = -0.25          # Lower point for cylinder expansion
rU =  0.25          # Upper point for cylinder expansion
sU =  rU + 400 * te   # Upper point for slot expansion
gU =  sU + 4 * tu   # Upper point for top plate expansion
delta_x = 1e-8

# Dependent Dimensions

wt = tu * 10             # Plenum wall thickness
hi = 2*Rc + te + wt     # Plenum height
le = .20 - Rc - hi/2    # Plenum exit x coordinate


 # Additional Coanda Dimensions:
lc = .02
ucx = le - lc
ucy = hi/2

phi = np.arctan((wt-tu)/lc) * 180 / np.pi



# print(fil,p0)

pback = np.zeros([22,2])

# Define Point Matrix:

pback[ 0,:] = np.array([le, Rc])                        # 0
pback[ 1,:] = np.array([le+delta_x, Rc + te])           # 2
pback[ 2,:] = np.array([le+2*delta_x, Rc + te + tu])    # 4
pback[ 3,:] = np.array([ucx, ucy])                      # 6
pback[ 4,:] = np.array([hi/2, hi/2])                    # 8
pback[ 5,:] = np.array([hi/2,-hi/2])                    # 10
pback[ 6,:] = np.array([ucx, -ucy])                     # 12
pback[ 7,:] = np.array([le+2*delta_x, -(Rc + te + tu)]) # 14
pback[ 8,:] = np.array([le+delta_x,-(Rc + te)])         # 16
pback[ 9,:] = np.array([le, -Rc])                       # 18
pback[10,:] = np.array([hi/2, bD])                      # 20
pback[11,:] = np.array([ucx, bD])                       # 22
pback[12,:] = np.array([le+3*delta_x,bD])               # 24
pback[13,:] = np.array([le+3*delta_x,bU])               # 26
pback[14,:] = np.array([ucx, bU])                       # 28
pback[15,:] = np.array([hi/2, bU])                      # 30
pback[16,:] = np.array([0, 0])                          # 32
pback[17,:] = np.array([-1, 0])                         # 34

# Finally: Define the blocking and grading parameters!

blocks_x_slot = 15
blocks_x_cont = 25
blocks_y_slot = 15
blocks_x_out = 46
blocks_y_out = 46
blocks_x_cyl = 75
blocks_y_cyl = 10

grade_x_flat = 1
grade_y_flat = 3000

grade_x_cont = 4

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
    '\t\thex (16  2  0 18 17  3  1 19) ({} {} 1) simpleGrading ( 1   1  1) // 6\n'.format(blocks_x_cyl, blocks_y_cyl),
    '\t\thex (14  4  2 16 15  5  3 17) ({} {} 1) simpleGrading ( 1   1  1) // 6\n'.format(blocks_x_cyl, blocks_y_cyl),
    '\t\thex (24 26  4 14 25 27  5 15) ({} {} 1) simpleGrading ( 1  {}  1) // 6\n'.format(blocks_x_cyl, blocks_y_out, 1/grade_y_flat),
    '\t\thex (26 28  6  4 27 29  7  5) ({} {} 1) simpleGrading ({}  {}  1) // 6\n'.format(blocks_x_cont, blocks_y_out, grade_x_cont, 1/grade_y_flat),
    '\t\thex (28 30  8  6 29 31  9  7) ({} {} 1) simpleGrading ( 1  {}  1) // 7\n'.format(blocks_x_out+14, blocks_y_out, 1/grade_y_flat),
    '\t\thex (30 34 32  8 31 35 33  9) ({} {} 1) simpleGrading ( 1  {}  1) // 8\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    '\t\thex (34 20 10 32 35 21 11 33) ({} {} 1) simpleGrading ( 1  {}  1) // 8\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    '\t\thex (20 22 12 10 21 23 13 11) ({} {} 1) simpleGrading ( 1  {}  1) // 9\n'.format(blocks_x_out+14, blocks_y_out, 1/grade_y_flat),
    '\t\thex (22 24 14 12 23 25 15 13) ({} {} 1) simpleGrading ({}  {}  1) //10\n'.format(blocks_x_cont, blocks_y_out, 1/grade_x_cont, 1/grade_y_flat),
    ');\n\n'
]

edges = [
    'edges\n(\n',
    '\t\tarc   0 18 180.0 (0 0 -1)\n',
    '\t\tarc   1 19 180.0 (0 0 -1)\n',
    '\t\tarc   2 16 180.0 (0 0 -1)\n',
    '\t\tarc   3 17 180.0 (0 0 -1)\n',
    '\t\tarc   4 14 180.0 (0 0 -1)\n',
    '\t\tarc   5 15 180.0 (0 0 -1)\n',
    '\t\tarc  30 34  90.0 (0 0  1)\n',
    '\t\tarc  31 35  90.0 (0 0  1)\n',
    '\t\tarc  34 20  90.0 (0 0  1)\n',
    '\t\tarc  35 21  90.0 (0 0  1)\n',
    '\t\tarc  32  8  90.0 (0 0 -1)\n',
    '\t\tarc  33  9  90.0 (0 0 -1)\n',
    '\t\tarc  32 10  90.0 (0 0  1)\n',
    '\t\tarc  33 11  90.0 (0 0  1)\n',
    '\t\tarc   4  6    {} (0 0  1)\n'.format(phi),
    '\t\tarc   5  7    {} (0 0  1)\n'.format(phi),
    '\t\tarc  12 14    {} (0 0  1)\n'.format(phi),
    '\t\tarc  13 15    {} (0 0  1)\n'.format(phi),
    '\t\tarc  26 24 180.0 (0 0 -1)\n',
    '\t\tarc  27 25 180.0 (0 0 -1)\n',
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
