# - readPyTreeFromPaths (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'test.adf')

t = Filter.convertFile2SkeletonTree('test.adf', maxDepth=3)

# Complete t par leur paths
Filter._readPyTreeFromPaths(t, 'test.adf', ['/Base/cart/GridCoordinates', 'Base/cart.0/GridCoordinates'])
test.testT(t, 1)
