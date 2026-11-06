# - copyValue (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])

# t2 has numpy copy only for nodes of type 'DataArray_t'
t2 = Internal.copyValue(t, byType='DataArray_t')

# t3 has numpy copy only for nodes of name 'Coordinate*'
t3 = Internal.copyValue(t, byName='Coordinate*')
