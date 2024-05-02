# - getByType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base',a])

# Return a standard node containing nodes of type 'Zone_t' as children
zones = Internal.getByType(t, 'Zone_t'); print(zones)
#>> ['Zone_t', None, [['cart', array([[10,  9,  0],...