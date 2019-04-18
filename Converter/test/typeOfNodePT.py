# - typeOfNode (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base1', 'Base2'])
print(Internal.typeOfNode(t))
#>> 3
