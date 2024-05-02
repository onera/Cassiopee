# - setName (pyTree) -
import Converter.Internal as Internal

node = Internal.createNode('node1', 'DataArray_t', value=1.)
Internal.setName(node, 'myNode'); print(node)
#>> ['myNode', array([ 1.]), [], 'DataArray_t']
