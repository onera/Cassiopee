# - addGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test
#
# 2D Case
#
a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 1))
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
a = C.initVars(a, '{F}={CoordinateY}*{CoordinateX}')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# partiellement coincident
a2 = G.cart((1., 0.4, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
a2 = T.oneovern(a2,(2,2,1))
a2 = C.initVars(a2, '{centers:Density}={centers:CoordinateX}')
a2 = C.initVars(a2, '{F}={CoordinateY}*{CoordinateX}')
a2 = C.addBC2Zone(a2, 'overlap1', 'BCOverlap', 'imax')
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, a2]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectNearMatch(t, dim=2)
t = Internal.addGhostCells(t,t,2,adaptBCs=0,fillCorner=1)
t = Internal.rmGhostCells(t,t, 2, adaptBCs=0)
test.testT(t,1)
#
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, a2]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectNearMatch(t, dim=2)
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=1)
t = Internal.rmGhostCells(t,t, 2, adaptBCs=1)
test.testT(t,2)

# geometrical extrapolation of corner cells
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, a2]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectNearMatch(t, dim=2)
t = Internal.addGhostCells(t,t,2,adaptBCs=0,fillCorner=0)
t = Internal.rmGhostCells(t,t,2,adaptBCs=0)
test.testT(t,3)

# geometrical extrapolation of corner cells
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, a2]
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = X.connectNearMatch(t, dim=2)
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=0)
t = Internal.rmGhostCells(t,t,2,adaptBCs=1)
test.testT(t,4)
