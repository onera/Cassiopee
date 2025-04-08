# - groupBCByType (PyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))

a = C.fillEmptyBCWith(a,'wall','BCWall')

t = C.newPyTree(['Base',a])
Internal._groupBCByBCType(t,'BCWall','FamWall')

test.testT(t,1)
