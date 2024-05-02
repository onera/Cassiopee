# - splitSize (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

# Structure + Champs + CL
a = G.cart((0,0,0),(1,1,1),(50,20,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
C._addVars(a, 'Density'); C._initVars(a, 'centers:cellN',1)
t = C.newPyTree(['Base', a])
t = T.splitSize(t, 1000)
t = X.connectMatch(t)
test.testT(t, 1)

# 2D Structure + Champs + CL
a = G.cart((0,0,0),(1,1,1),(50,20,2))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
C._addVars(a, 'Density'); C._initVars(a, 'centers:cellN',1)
t = C.newPyTree(['Base',3,a])
t = T.splitSize(t, 100)
t = X.connectMatch(t)
test.testT(t,2)
