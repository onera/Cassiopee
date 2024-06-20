# - convertPyTree2File (pyTree) -
# - openfoam -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# NGON3
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
C._fillEmptyBCWith(z, 'far', 'BCFarfield')
C.convertPyTree2File(z, LOCAL+'/out1.foam')
test.testF(LOCAL+'/out1.foam/constant/polyMesh/faces', 1)
test.testF(LOCAL+'/out1.foam/constant/polyMesh/owner', 2)
test.testF(LOCAL+'/out1.foam/constant/polyMesh/points', 3)
test.testF(LOCAL+'/out1.foam/constant/polyMesh/neighbour', 4)

# NGON4
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
C._fillEmptyBCWith(z, 'far', 'BCFarfield')
C.convertPyTree2File(z, LOCAL+'/out2.foam')
test.testF(LOCAL+'/out2.foam/constant/polyMesh/faces', 1)
test.testF(LOCAL+'/out2.foam/constant/polyMesh/owner', 2)
test.testF(LOCAL+'/out2.foam/constant/polyMesh/points', 3)
test.testF(LOCAL+'/out2.foam/constant/polyMesh/neighbour', 4)
