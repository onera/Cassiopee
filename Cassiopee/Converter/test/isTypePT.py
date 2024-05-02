# - isType (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base1', 'Base2'])

# This is true
print(Internal.isType(t, 'CGNSTree_t'))
#>> True

# This is false
print(Internal.isType(t, 'CGNSBase_t'))
#>> False

# Check with wildcard: answer true if the type matches the string
print(Internal.isType(t, 'CGNS*'))
#>> True
