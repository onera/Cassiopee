# - getNodePosition (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
g = Internal.getNodeFromType(a, 'GridCoordinates_t')
x = Internal.getNodeFromName(g, 'CoordinateZ')

print(Internal.getNodePosition(x, g))
#>> 2
