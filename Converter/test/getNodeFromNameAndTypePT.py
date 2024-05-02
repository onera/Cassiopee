# - getNodeFromNameAndType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# Return the node named 'cart' and type 'Zone_t'
node = Internal.getNodeFromNameAndType(t, 'cart', 'Zone_t'); print(node)
#>> ['cart', array([..]), [..], 'Zone_t']
