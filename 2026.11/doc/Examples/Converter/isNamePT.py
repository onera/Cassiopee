# - isName (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base1', 'Base2'])

# This is false
print(Internal.isName(t[2][1], 'Base3'))
#>> False

# This is true if node name matches the string
print(Internal.isName(t[2][1], 'Base*'))
#>> True
