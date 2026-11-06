# - readNodesFromPaths (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, LOCAL+'/test.hdf')

# Relit les noeuds par leur paths
nodes = Filter.readNodesFromPaths(LOCAL+'/test.hdf', '/Base/cart/GridCoordinates')
test.testT(nodes, 1)
