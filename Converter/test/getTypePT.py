# - getType (pyTree) -
import Converter.Internal as Internal

node = Internal.createNode('myNode', 'DataArray_t', value=1.)
print(Internal.getType(node))
#>> DataArray_t
