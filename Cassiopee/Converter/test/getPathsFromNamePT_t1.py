# - getPathsFromName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# Return paths of zone named 'cart'
paths = Internal.getPathsFromName(t, 'cart')
test.testO(paths, 1)

# Return the 3 coordinate paths
paths = Internal.getPathsFromName(t, 'Coordinate*')
test.testO(paths, 2)
