# - getNodePosition (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
g = Internal.getNodeFromType1(a, 'GridCoordinates_t')
x = Internal.getNodeFromName1(g, 'CoordinateZ')

test.testO(Internal.getNodePosition(x, g))
