# - getNodesFromValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base',a])

# Return nodes with given value
nodes = Internal.getNodesFromValue(t, 2.4); print(nodes)
#>> []
nodes = Internal.getNodesFromValue(t, 'Structured'); print(nodes)
#>> [['ZoneType', array(..), [], 'ZoneType_t']]
