# - getPathsFromName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# Return paths of zone named 'cart'
paths = Internal.getPathsFromName(t, 'cart'); print(paths)
#>> ['CGNSTree/Base/cart']

# Return the 3 coordinate paths
paths = Internal.getPathsFromName(t, 'Coordinate*'); print(paths)
#>> ['CGNSTree/Base/cart/GridCoordinates/CoordinateX', '/CGNSTree/Base/cart/GridCoordinates/CoordinateY', '/CGNSTree/Base/cart/GridCoordinates/CoordinateZ']
