# - checkPyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cart( (9,0,0), (1,1,1), (10,10,10) )
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = X.connectMatch(t)

# check nodes conformity
errors = Internal.checkPyTree(t, level=1)
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
# check types
errors += Internal.checkPyTree(t, level=8)
test.testO(errors, 1)
