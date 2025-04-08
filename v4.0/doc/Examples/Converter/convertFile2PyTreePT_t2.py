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
a = C.convertFile2PyTree(LOCAL+'/out1.foam')
test.testT(a, 1)

# NGON4
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
C._fillEmptyBCWith(z, 'far', 'BCFarfield')
C.convertPyTree2File(z, LOCAL+'/out2.foam')
a = C.convertFile2PyTree(LOCAL+'/out2.foam')
test.testT(a, 2)
