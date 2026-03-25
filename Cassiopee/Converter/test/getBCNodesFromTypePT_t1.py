# - getBCNodesFromType (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (10,10,2))
C._addBC2Zone(a, 'wall1', 'BCWallInviscid', 'jmin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
C._addBC2Zone(a, 'wall2', 'BCWallViscous', 'jmax')
C._addBC2Zone(a, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(a, 'outlet', 'BCOutflow', 'kmax')
t = C.newPyTree(['Base', a])
bc = Internal.getNodeFromType2(a, 'BC_t')

# On a zone
n = Internal.getBCNodesFromType(a)  # all BCs
test.testO(n, 1)

n = Internal.getBCNodesFromType(a, 'BCSymmetryPlane')  # not found, empty list
test.testO(n, 2)

n = Internal.getBCNodesFromType(a, 'BCWallInviscid')  # one BC
test.testO(n, 3)

n = Internal.getBCNodesFromType(a, ['BCInflow', 'BCOutflow'])  # from a list of BC types
test.testO(n, 4)

n = Internal.getBCNodesFromType(a, 'BCWall*')  # wildcard
test.testO(n, 5)

# On a PyTree
n = Internal.getBCNodesFromType(t)  # all BCs
test.testO(n, 6)

n = Internal.getBCNodesFromType(t, 'BCSymmetryPlane')  # not found, empty list
test.testO(n, 7)

n = Internal.getBCNodesFromType(t, 'BCWallInviscid')  # one BC
test.testO(n, 8)

n = Internal.getBCNodesFromType(t, 'BCWall*')  # wildcard
test.testO(n, 9)

# On a BC_t node
n = Internal.getBCNodesFromType(bc)  # all BCs
test.testO(n, 10)

n = Internal.getBCNodesFromType(bc, 'BCSymmetryPlane')  # not found, empty list
test.testO(n, 11)

n = Internal.getBCNodesFromType(bc, 'BCWallInviscid')  # one BC
test.testO(n, 12)

n = Internal.getBCNodesFromType(bc, 'BCWall*')  # wildcard
test.testO(n, 13)

# On a list of zones
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((5,0.9,0), (0.01,0.01,1.), (20,20,2)); b[0] = 'cart2'
zones = [a, b]
C._addBC2Zone(zones, 'wall1', 'BCWallInviscid', 'imin')
C._addBC2Zone(zones, 'wall2', 'BCWallViscous', 'jmin')
C._addBC2Zone(zones, 'wall3', 'BCWallViscousIsothermal', 'jmax')
C._addBC2Zone(zones, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(zones, 'farfield', 'BCFarfield', 'kmax')

n = Internal.getBCNodesFromType(zones)  # all BCs
test.testO(n, 14)

n = Internal.getBCNodesFromType(zones, 'BCSymmetryPlane')  # not found, empty list
test.testO(n, 15)

n = Internal.getBCNodesFromType(zones, 'BCWallViscous')  # one BC
test.testO(n, 16)

n = Internal.getBCNodesFromType(zones, 'BCWallViscous*')  # wildcard
test.testO(n, 17)
