# - addGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test
#
# 3D Case
y0 = 0.; dimPb = 3; nk = 3
a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21,nk))
C._initVars(a, '{centers:Density}={centers:CoordinateX}')
C._initVars(a, '{F}={CoordinateY}*{CoordinateX}')
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# partiellement coincident
a2 = G.cart((1.,y0, 0.), (0.1, 0.1, 0.1), (11, 21,nk))
a2 = T.oneovern(a2,(2,2,1))
C._initVars(a2, '{centers:Density}={centers:CoordinateX}')
C._initVars(a2, '{F}={CoordinateY}*{CoordinateX}')
C._addBC2Zone(a2, 'overlap1', 'BCOverlap', 'imax')
t = C.newPyTree(['Base',a,a2])
t = X.connectNearMatch(t,2,dim=dimPb)
t = Internal.addGhostCells(t,t,2,adaptBCs=0,fillCorner=1)
test.testT(t,1)
#
t = C.newPyTree(['Base',a,a2])
C._addState(t[2][1], 'EquationDimension',dimPb)
t = X.connectNearMatch(t,2)
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=1)
test.testT(t,2)

# geometrical extrapolation of corner cells
t = C.newPyTree(['Base',a,a2])
C._addState(t[2][1], 'EquationDimension',dimPb)
t = X.connectNearMatch(t,2)
t = Internal.addGhostCells(t,t,2,adaptBCs=0,fillCorner=0)
test.testT(t,3)

# geometrical extrapolation of corner cells
t = C.newPyTree(['Base',a,a2])
C._addState(t[2][1], 'EquationDimension',dimPb)
t = X.connectNearMatch(t,2)
t = Internal.addGhostCells(t,t,2,adaptBCs=1,fillCorner=0)
test.testT(t,4)
