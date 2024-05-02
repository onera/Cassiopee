# - isTopTree (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base1', 'Base2'])
print(Internal.isTopTree(t))
#>> True
print(Internal.isTopTree(t[2][1]))
#>> False
