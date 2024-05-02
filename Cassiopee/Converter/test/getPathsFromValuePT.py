# - getPathsFromValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# Return nodes with given value
paths = Internal.getPathsFromValue(t, 3.1); print(paths)
#>> ['CGNSTree/CGNSLibraryVersion']
paths = Internal.getPathsFromValue(t, 'Structured'); print(paths)
#>> ['CGNSTree/Base/cart/ZoneType']
