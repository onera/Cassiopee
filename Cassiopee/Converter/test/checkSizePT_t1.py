# - checkSize (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test
a = G.cart((0,0,0), (1,1,1), (11,11,11))
b = G.cart((10,0,0), (1,1,1), (11,11,11))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base',a,b])
t = X.connectMatch(t)

# check nodes conformity
errors = Internal.checkSize(t, sizeMax=500)
test.testO(errors, 1)
