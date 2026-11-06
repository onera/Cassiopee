# - getNobOfBase (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base', 'Base2'])
b = Internal.getNodeFromName(t, 'Base2')
nob = C.getNobOfBase(b, t); print(nob)
#>> 2
# This means that t[2][nob] = b
