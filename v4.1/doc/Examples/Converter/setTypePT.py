# - setType (pyTree) -
import Converter.Internal as Internal

node = Internal.createNode('node1', 'DataArray_t', value=1.)
Internal.setType(node, 'Zone_t'); print(node)
#>> ['node1', array([ 1.]), [], 'Zone_t']
