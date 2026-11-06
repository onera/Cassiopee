# - copyValue (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (5,5,5))
t = C.newPyTree(['Base']); t[2][1][2] += [a]

# t2 has numpy copy only for nodes of type 'DataArray_t'
t2 = Internal.copyValue(t, byType='DataArray_t')
x = Internal.getNodeFromName(t2, 'CoordinateX')
x[1][0,0,0] = 5.
test.testT(t2, 1)
