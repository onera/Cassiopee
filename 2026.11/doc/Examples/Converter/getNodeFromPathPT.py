# - getNodeFromPath (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])

# Return GridCoordinates node
coords = Internal.getNodeFromPath(t, 'CGNSTree/Base/cart/GridCoordinates'); print(coords)
#>> ['GridCoordinates', None, [..], 'DataArray_t']

# Return GridCoordinates node (path is relative to input node)
coords = Internal.getNodeFromPath(a, 'cart/GridCoordinates'); print(coords)
#>> ['GridCoordinates', None, [..], 'DataArray_t']

# Return GridCoordinates node (path is relative to input node)
coords = Internal.getNodeFromPath(a, './GridCoordinates'); print(coords)
#>> ['GridCoordinates', None, [..], 'DataArray_t']
