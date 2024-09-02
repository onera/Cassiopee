# - isChild (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base1', 'Base2'])

b1 = Internal.getNodeFromPath(t, 'CGNSTree/Base1')
b2 = Internal.getNodeFromPath(t, 'CGNSTree/Base2')
n = Internal.createChild(b1, 'mynode', 'DataArray_t', value=1.)

print(Internal.isChild(b1, n))
#>> True
print(Internal.isChild(b2, n))
#>> False
