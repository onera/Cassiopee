# - deleteEmptyZones (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (3,3,3))
b = P.selectCells(a, '{CoordinateX} > 12')

# Sur un pyTree
tp = C.newPyTree(['Base',a,b])
tp = C.deleteEmptyZones(tp)
test.testT(tp, 1)

# Sur une liste de zones
l = C.deleteEmptyZones([a,b])
test.testT(l, 2)
