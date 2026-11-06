# - getNodesFromNameAndType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# Return nodes named 'cart' and of type 'Zone_t'
nodes = Internal.getNodesFromNameAndType(t, 'cart', 'Zone_t'); print(nodes)
#>> [['cart', array([..]), [..], 'Zone_t']]

# Return the 3 coordinate nodes
nodes = Internal.getNodesFromNameAndType(t, 'Coordinate*', 'Data*_t'); print(nodes)
#>> [['CoordinateX', array([..]), [], 'DataArray_t'], ['CoordinateY', array([]), 'DataArray_t'], ['CoordinateX', array([..]), [], 'DataArray_t']]
