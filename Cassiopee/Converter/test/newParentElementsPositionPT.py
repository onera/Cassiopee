# - newParentElementsPosition (pyTree) -
import Converter.Internal as Internal

# Create a parent elements position node
n = Internal.newParentElementsPosition(value=None); Internal.printTree(n)
#>> ['ParentElementsPosition',None,[0 son],'DataArray_t']

# Attach it to elements
d = Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, eboundary=0)
Internal.newParentElementsPosition(None, parent=d)
