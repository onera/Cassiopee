# - getNodeFromType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base',a])

# Return the first node of type 'Zone_t'
node = Internal.getNodeFromType(t, 'Zone_t'); print(node)

# Limit search to second level (faster)
node = Internal.getNodeFromType2(t, 'Zone_t'); print(node)
#>> ['cart', array([[10,  9,  0],...