# - deletePaths (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])
C.convertPyTree2File(t, LOCAL+'/out.hdf')

# Delete paths
Filter.deletePaths(LOCAL+'/out.hdf', 'CGNSTree/Base/cart/GridCoordinates')

# Reread for test
t = C.convertFile2PyTree(LOCAL+'/out.hdf')
test.testT(t, 1)
