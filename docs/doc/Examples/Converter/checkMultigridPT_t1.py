# - checkMultigrid (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (11,11,10) )
b = G.cart( (10,0,0), (1,1,1), (11,11,10) )
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = X.connectMatch(t)

# check multigrid compatibility
errors = Internal.checkMultigrid(t, level=1, nbMinCoarseW=3, nbMinCoarseB=5)
test.testO(errors, 1)

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
b = G.cart( (10,0,0), (1,1,1), (11,11,11) )
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = X.connectMatch(t)

# check nodes conformity
errors = Internal.checkMultigrid(t, level=1, nbMinCoarseW=3, nbMinCoarseB=5)
test.testO(errors,2)

a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
b = G.cart( (10,8,0), (1,1,1), (11,11,11) )
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = X.connectMatch(t)
# check nodes conformity
errors = Internal.checkMultigrid(t, level=1, nbMinCoarseW=3, nbMinCoarseB=5)
test.testO(errors, 3)
