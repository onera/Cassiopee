# - checkPyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((9,0,0), (1,1,1), (10,10,10))
c = G.cartTetra((9,9,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base',a,b,c])
t = X.connectMatch(t)

errors = []
# check unique base names
errors += Internal.checkPyTree(t, level=2)
# check unique zone names
errors += Internal.checkPyTree(t, level=3)
# check unique BC names
errors += Internal.checkPyTree(t, level=4)
# check BC ranges
errors += Internal.checkPyTree(t, level=5)
# check opposite ranges
errors += Internal.checkPyTree(t, level=6)
# check family definition
errors += Internal.checkPyTree(t, level=7)
# check CGNSTypes
errors += Internal.checkPyTree(t, level=8)
# check element nodes
errors += Internal.checkPyTree(t, level=9); print(errors)
#>> []

# Introduce errors (on purpose!)
n = Internal.getNodeFromType(t, 'Zone_t')
Internal.setType(n, 'Zon_t')
errors = Internal.checkPyTree(t, level=8); print(errors)
#>> [['cart', array(..), [..], 'Zon_t'], 'Unknown CGNS type Zon_t for node cart.\n']
