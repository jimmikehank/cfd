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


def store(retain,name):
    import os
    import shutil
    ignore = ['data_process.ipynb','autoMesh.py','.ipynb_checkpoints']
    copy = ['0','system','constant','dynamicCode']
    dirs = os.listdir()
    delete = []
    for item in dirs:
        if item[0:3] == 'pro':
            delete.append(item)
    target = name
    # casefile = "/home/james/Documents/research/completed_cases/coanda_airfoils/{}/".format(target)
    casefile = "/media/james/Data/james/completed_cases/coanda_plenum_comparison/{}/".format(target)
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

# This function parses input arguments to generate mesh allowing external programs to iterate through different mesh configurations

def arg_handle(args):
    # Default Values
    Rc = 0.003556   # Coanda cylinder radius
    te = 0.0001778  # R/te = 20
    tu = 0.00025  # Upper surface thickness at exit

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
        elif current == '-save':
            store(retain,args[i+1])
    return Rc, te, tu

Rc, te, tu = arg_handle(args)



#-------------------------------------------------------------#

print('\nCoanda Surface Radius = {}\nSlot Exit Width = {}\nExit Surface Thickness = {}\n'.format(Rc, te, tu))

# blockMeshDefaults
scale = 1


# Bounding Box:
farfield = 5
bU =  farfield # Upper bound
bD = -farfield # Lower bound
bL = -farfield # Left bound
bR =  farfield # Right bound
fr =  0.078 # Front bound
bk = -0.078 # Back bound

# Mesh expansion for coanda cylinder

rL = -0.25          # Lower point for cylinder expansion
rU =  0.25          # Upper point for cylinder expansion
sU =  rU + 400 * te   # Upper point for slot expansion
gU =  sU + 4 * tu   # Upper point for top plate expansion
delta_x = 1e-8

# Dependent Dimensions

wt = tu * 10             # Plenum wall thickness
hi = 2*Rc + te + wt     # Plenum height
le = .3 - Rc -hi/2  # Plenum exit x coordinate


 # Additional Coanda Dimensions:
lc = .02
ucx = le - lc
ucy = hi/2

phi = np.arctan((wt-tu)/lc) * 180 / np.pi

# New coordinates for Plenum:
pln = 0.04              # Plenum length
pcn = (Rc - te)*1.5
hcn = (Rc-te)/2
beta = np.arctan(pcn/hcn) * 180 / np.pi
ai = 180 - 2*beta
# print(fil,p0)

pback = np.zeros([30,2])

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
pback[17,:] = np.array([bL, 0])                         # 34

# New points to define Plenum
pback[18,:] = np.array([-pln+le, Rc+te])
pback[19,:] = np.array([-pln+le, te])
pback[20,:] = np.array([-2*pcn+le,te])
pback[21,:] = np.array([-2*pcn+le,Rc+te])
pback[22,:] = np.array([-pcn+le,hcn+te])
pback[23,:] = np.array([-pcn+le,Rc+te])

pback[24,:] = np.array([-pln+le, -(Rc+te)])
pback[25,:] = np.array([-pln+le, -te])
pback[26,:] = np.array([-2*pcn+le,-te])
pback[27,:] = np.array([-2*pcn+le,-(Rc+te)])
pback[28,:] = np.array([-pcn+le,-(hcn+te)])
pback[29,:] = np.array([-pcn+le,-(Rc+te)])

# Finally: Define the blocking and grading parameters!

blocks_x_slot = 15
blocks_x_cont = 25
blocks_y_slot = 10
blocks_x_out = 70
blocks_y_out = 80
blocks_x_cyl = 75
blocks_y_cyl = 15

blocks_x_in = 20
blocks_y_in = blocks_y_cyl

grade_x_flat = 1
grade_y_flat = 2000

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
    '\t\thex (28 30  8  6 29 31  9  7) ({} {} 1) simpleGrading (((0.5 0.5  2) (0.5 0.5 0.1))  {}  1) // 7\n'.format(blocks_x_out+50, blocks_y_out, 1/grade_y_flat),
    '\t\thex (30 34 32  8 31 35 33  9) ({} {} 1) simpleGrading ( 1  {}  1) // 8\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    '\t\thex (34 20 10 32 35 21 11 33) ({} {} 1) simpleGrading ( 1  {}  1) // 8\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    '\t\thex (20 22 12 10 21 23 13 11) ({} {} 1) simpleGrading (((0.5 0.5 10) (0.5 0.5 0.5))  {}  1) // 9\n'.format(blocks_x_out+50, blocks_y_out, 1/grade_y_flat),
    '\t\thex (22 24 14 12 23 25 15 13) ({} {} 1) simpleGrading ({}  {}  1) //10\n'.format(blocks_x_cont, blocks_y_out, 1/grade_x_cont, 1/grade_y_flat),
    '\t\thex (38 40 42 36 39 41 43 37) ({} {} 1) simplegrading ( 1   1  1) //11\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (40 44 46 42 41 45 47 43) ({} {} 1) simpleGrading ( 1   1  1) //12\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (44  0  2 46 45  1  3 47) ({} {} 1) simpleGrading ( 1   1  1) //13\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (48 54 52 50 49 55 53 51) ({} {} 1) simpleGrading ( 1   1  1) //14\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (54 58 56 52 55 59 57 53) ({} {} 1) simpleGrading ( 1   1  1) //15\n'.format(blocks_x_in,blocks_y_in),
    '\t\thex (58 16 18 56 59 17 19 57) ({} {} 1) simpleGrading ( 1   1  1) //16\n'.format(blocks_x_in,blocks_y_in),
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
    '\t\tarc  40 44    {} (0 0  1)\n'.format(ai),
    '\t\tarc  44  0    {} (0 0 -1)\n'.format(ai),
    '\t\tarc  41 45    {} (0 0  1)\n'.format(ai),
    '\t\tarc  45  1    {} (0 0 -1)\n'.format(ai),

    '\t\tarc  52 56    {} (0 0 -1)\n'.format(ai),
    '\t\tarc  56 18    {} (0 0  1)\n'.format(ai),
    '\t\tarc  53 57    {} (0 0 -1)\n'.format(ai),
    '\t\tarc  57 19    {} (0 0  1)\n'.format(ai),
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