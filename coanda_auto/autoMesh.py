import numpy as np
import matplotlib.pyplot as plt

# blockMeshDefaults
scale = 1

# Control Variables:

Rc = 0.01   # Coanda cylinder radius
te = 0.001  # R/te = 10
tu = 0.002  # Upper surface thickness at exit
ru = 0.50    # Contour radius of upper surface to exit
ai = 60.0   # Contraction angle inside plenum

# Bounding Box:

bU =  1.0 # Upper bound
bD = -1.0 # Lower bound
bL = -1.0     # Left bound
bR =  1.0 # Right bound
fr =  0.5   # Front bound
bk = -0.5   # Back bound

# Additional Coanda Dimensions:

pl =  0.10  # Plenum internal length
pi =  0.03  # Plenum additional mesh length
pw =  0.008 # Plenum wall thickness

# Mesh expansion for coanda cylinder

rL = -0.1            # Lower point for cylinder expansion
rU =  0.1            # Upper point for cylinder expansion
sU =  rU + 4 * te    # Upper point for slot expansion
gU =  sU + 2 * tu    # Upper point for top plate expansion

# Dependent Dimensions

hi = 2*Rc + te - pw  # Plenum height
us = hi + pw         # Upper surface of plenum
eu = hi + tu         # Upper surface of slot exit
le = pl + Rc         # Plenum exit x coordinate
ly = hi - te         # Plenum exit y coordinate
pm = (hi)/2          # Middle of plenum for leading edge curve
pn = -(hi/2 + pw)    # LE of curved surface
CC = Rc - pw         # Cylinder center height
ri = (CC + Rc * np.cos(np.radians(ai)))/(1 - np.cos(np.radians(ai)))

# Fillet calculations:

fil = np.zeros([2,2])

fil[0,0] = le - (np.sin(np.radians(ai)) * (Rc+ri))
fil[0,1] = 0
fil[1,0] = le - Rc * np.sin(np.radians(ai))
fil[1,1] = CC + Rc * np.cos(np.radians(ai))
#print(fil)
topang = np.degrees(np.arctan(abs(us-eu)/abs(le-fil[0,0]))) # Fillet for upper surface to exit

pback = np.zeros([26,2])

# Define Point Matrix:

pback[ 0,:] = np.array([0,0])
pback[ 1,:] = fil[0,:]
pback[ 2,:] = np.array([pi,hi])
pback[ 3,:] = np.array([0,hi])
pback[ 4,:] = fil[1,:]
pback[ 5,:] = np.array([pl,hi])
pback[ 6,:] = np.array([le,ly])
pback[ 7,:] = np.array([le,hi])
pback[ 8,:] = np.array([le,-pw])
pback[ 9,:] = np.array([0,-pw])
pback[10,:] = np.array([0,bD])
pback[11,:] = np.array([le,bD])
pback[12,:] = np.array([bR,bD])
pback[13,:] = np.array([bR,rL])
pback[14,:] = np.array([bR,rU])
pback[15,:] = np.array([bR,sU])
pback[16,:] = np.array([bR,gU])
pback[17,:] = np.array([bR,bU])
pback[18,:] = np.array([le,bU])
pback[19,:] = np.array([fil[0,0],bU])
pback[20,:] = np.array([0,bU])
pback[21,:] = np.array([0,us])
pback[22,:] = np.array([fil[0,0],us])
pback[23,:] = np.array([le,eu])
pback[24,:] = np.array([bL,pm])
pback[25,:] = np.array([pn,pm])

print(pback[0:10,:])
plt.plot(pback[0:10,0],pback[0:10,1])
plt.show()
#pback = np.around(pback,5)


header = ['/*---------------------------------*- C++ -*-----------------------------------*/\n','//    //////////////      ///////////////\n','//    //////////////      //////////////\n','//         ////           ////             /////       ////  //  /// //// ////\n','//         ////     ////  ///////////   ///    ///  ///    ///   ////  ///  ///\n','//  ////   ////           ///////////   ///    ///  ///    ///   ///   ///  ///\n','//   //// ////            ////          ///    ///  ///    ///   ///   ///  ///\n','//     /////              ////            /////       ///// ///  ///   ///  ///\n','/*---------------------------------*- C++ -*-----------------------------------*/\n\n']

opener = ['FoamFile','\n','{','\n','version\t\t\t2.0;','\n','format\t\t\tascii;','\n','class\t\t\t\tdictionary;','\n','object\t\t\tblockMeshDict;','\n','}\n\n','// ***************************************************************************** //\n\n']

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

blocks = ['blocks\n(\n','\t\thex ( 0  2  4  6  1  3  5  7) ( 50  80 1) simpleGrading (  1     1  1) // 0\n','\t\thex ( 2  8 10  4  3  9 11  5) ( 40  80 1) simpleGrading (.25     1  1) // 1\n','\t\thex ( 8 12 14 10  9 13 15 11) (120  80 1) simpleGrading (  1     1  1) // 2\n','\t\thex (20 22 16 18 21 23 17 19) (300  50 1) simpleGrading ( .5  .004  1) // 3\n','\t\thex (24 26 16 22 25 27 17 23) ( 50 100 1) simpleGrading (.004 .004  1) // 4\n','\t\thex (26 28 12 16 27 29 13 17) (150 100 1) simpleGrading (  1  .004  1) // 5\n','\t\thex (28 30 14 12 29 31 15 13) ( 80 100 1) simpleGrading (  1  .004  1) // 6\n','\t\thex (30 32 46 14 31 33 47 15) (100 100 1) simpleGrading (  1  .004  1) // 7\n','\t\thex (32 34 36 46 33 35 37 47) ( 50 100 1) simpleGrading (250  .004  1) // 8\n','\t\thex (36 38 44 46 37 39 45 47) ( 50  50 1) simpleGrading (  1  .004  1) // 9\n','\t\thex (38 40 42 44 39 41 43 45) ( 90  50 1) simpleGrading (  1  .004  1) // 10\n','\t\thex (40 48 50 42 41 49 51 43) (110  50 1) simpleGrading (  1  .004  1) // 11\n','\t\thex (48 20 18 50 49 21 19 51) (110  50 1) simpleGrading (  1  .004  1) // 12\n',');\n\n']

edges = ['edges\n(\n','\t\tarc   2  8  {} (0 0  1)\n'.format(ai),'\t\tarc   3  9  {} (0 0  1)\n'.format(ai),'\t\tarc   8 12  {} (0 0 -1)\n'.format(ai),'\t\tarc   9 13  {} (0 0 -1)\n'.format(ai),'\t\tarc  12 16 180.0 (0 0 -1)\n','\t\tarc  13 17 180.0 (0 0 -1)\n','\t\tarc  44 46 {} (0 0 -1)\n'.format(topang),'\t\tarc  45 47 {} (0 0 -1)\n'.format(topang),'\t\tarc  40 48  90.0 (0 0  1)\n','\t\tarc  41 49  90.0 (0 0  1)\n','\t\tarc  48 20  90.0 (0 0  1)\n','\t\tarc  49 21  90.0 (0 0  1)\n','\t\tarc  42 50  90.0 (0 0  1)\n','\t\tarc  43 51  90.0 (0 0  1)\n','\t\tarc  50 18  90.0 (0 0  1)\n','\t\tarc  51 19  90.0 (0 0  1)\n',');\n\n']

pointsdef = makepoints(pback,bk,fr)


with open('./system/blockMeshDict','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(scale)
    f.writelines(pointsdef)
    f.writelines(blocks)
    f.writelines(edges)
    with open('faces','r') as g:
        f.writelines(g.readlines())
