import Converter.PyTree as C
import Generator.PyTree as G
#import Post.PyTree as P
import os, sys
import Intersector.PyTree as XOR
import numpy


if len(sys.argv) is not 2 :
    print "ARG ERROR : 1 argument to provide : mesh to adapt"
    sys.exit()

ifile1=sys.argv[1] # mesh to adapt in NGON format
z = C.convertFile2PyTree(ifile1)

z = C.fillEmptyBCWith(z, 'wall', 'BCWall')


########################## create the hook
hmesh = XOR.createHMesh(z, 0) # 0 : ISOTROPIC subdivision 
######################################## 
    
for i in range(5): # simple loop to demonstrate the feature
  
  # nodal specification
  n = C.getNPts(z)
  nodal_vals = numpy.empty((n,), dtype=numpy.int32)
  nodal_vals[:] = 1
  # one nodal_vals and one hmesh per zone
  z = XOR.adaptCellsNodal(z, [nodal_vals], hmesh)

#C.convertPyTree2File(z, "hmesh.cgns") # all the hierarchy is in !

z = XOR.conformizeHMesh(z, hmesh)     # each children faces replace its mother in any polyhedron

z = XOR.closeOctalCells(z)            # close cells (adding point on lateral faces)

#C.convertPyTree2File(z, "out1.cgns")

########################## free the hook
XOR.deleteHMesh(hmesh);
#####################################


