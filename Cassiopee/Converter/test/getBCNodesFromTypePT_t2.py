# - getBCNodesFromType (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (10,10,2))
C._addBC2Zone(a, 'wall1', 'FamilySpecified:CYLINDER', 'jmin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
C._addBC2Zone(a, 'wall2', 'FamilySpecified:CYLINDER', 'jmax')
C._addBC2Zone(a, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(a, 'outlet', 'BCOutflow', 'kmax')
t = C.newPyTree(['Base', a])
C._addFamily2Base(t[2][1], 'CYLINDER', bndType='BCWallViscous')

n = Internal.getBCNodesFromType(t)  # all BCs
test.testO(n, 1)

n = Internal.getBCNodesFromType(t, 'BCWall*')  # wildcard
test.testO(n, 2)

n = Internal.getBCNodesFromType(t, 'BCWallViscous')  # by BC type
test.testO(n, 3)
