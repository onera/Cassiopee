# - getNodesFromType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])

# Return nodes of type 'Zone_t'
zones = Internal.getNodesFromType(t, 'Zone_t'); print(zones)
#>> [['cart', array(..), [..], 'Zone_t']]

# Limit search to 2 levels (faster)
zones = Internal.getNodesFromType2(t, 'Zone_t'); print(zones)
#>> [['cart', array(..), [..], 'Zone_t']]
