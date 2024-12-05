import Converter.PyTree as C
import Generator.PyTree as G
#import Post.PyTree as P
import os, sys
import Intersector.PyTree as XOR
import numpy


if len(sys.argv) is not 3:
    print("ARG ERROR : 2 arguments to provide : mesh source")
    sys.exit()

ifile1=sys.argv[1] # mesh to adapt in NGON format
ifile2=sys.argv[2] # source mesh
z = C.convertFile2PyTree(ifile1)
s = C.convertFile2PyTree(ifile2)

z = C.fillEmptyBCWith(z, 'wall', 'BCWall')


########################## create the hook
hmsh = XOR.createHMesh(z, 0) # 0 : ISOTROPIC subdivision
########################################

for i in range(5): # simple loop to demonstrate the feature

    #nodal specification
    n = C.getNPts(z)
    nodal_vals = numpy.empty((n,), dtype=Internal.E_NpyInt)
    nodal_vals[:] = 1
    #one nodal_vals and one hmesh per zone
    z = XOR.adaptCellsNodal(z, [nodal_vals], hmesh=hmsh)
    # refine now with source mesh
    z = XOR.adaptCells(z, s, itermax=-1, hmesh=hmsh)

#C.convertPyTree2File(z, "hmesh.cgns") # all the hierarchy is in !

z = XOR.conformizeHMesh(z, hmsh)     # each children faces replace its mother in any polyhedron

z = XOR.closeCells(z)            # close cells (adding point on lateral faces)

#C.convertPyTree2File(z, "out1.cgns")

########################## free the hook
XOR.deleteHMesh(hmsh);
#####################################
