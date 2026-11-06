# - correctPyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test

# Correct
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cart( (9,0,0), (1,1,1), (10,10,10) )
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = X.connectMatch(t)
print(Internal.checkPyTree(t))
t = Internal.correctPyTree(t)
test.testT(t, 1)

# Wrong range for BCWall and wrong range for donor
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cart( (0,0,0), (1,1,1), (10,10,10) ); b[0] = 'toto'
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
a = C.addBC2Zone(a, 'wall', 'BCWall', wrange=[1,1,1,11,1,10])
a = C.addBC2Zone(a, 'match', 'BCMatch', wrange=[1,1,1,10,1,10],
                 zoneDonor=b, rangeDonor=[1,1,1,10,1,11], trirac=[1,2,3])

t = C.newPyTree(['Base', a])
t = Internal.correctPyTree(t)
test.testT(t, 2)
