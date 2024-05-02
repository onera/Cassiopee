# - newParentElements (pyTree) -
import Converter.Internal as Internal

# Create a parent elements node
n = Internal.newParentElements(value=None); Internal.printTree(n)
#>> ['ParentElements',None,[0 son],'DataArray_t']

# Attach it to elements
d = Internal.newElements('Elements', 'UserDefined', None, 0)
Internal.newParentElements(None, parent=d)
