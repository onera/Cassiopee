# - getBCFaceNode (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# structure
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName(a, 'wall')
ind = Internal.getBCFaceNode(a, b)
test.testO(ind, 1)
