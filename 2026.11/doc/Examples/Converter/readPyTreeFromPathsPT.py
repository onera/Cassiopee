# - readPyTreeFromPaths (pyTree) -
import Converter.PyTree as C
import Converter.Filter as Filter
import Generator.PyTree as G

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'test.hdf')

t = Filter.convertFile2SkeletonTree('test.hdf', maxDepth=3)

# Complete t par leur paths
Filter._readPyTreeFromPaths(t, 'test.hdf', ['/Base/cart/GridCoordinates', 'Base/cart.0/GridCoordinates'])
C.convertPyTree2File(t, 'out.hdf')
