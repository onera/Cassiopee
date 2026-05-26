# - getBCNodesFromName (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (10,10,2))
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
C._addBC2Zone(a, 'wall2', 'BCWall', 'jmax')
C._addBC2Zone(a, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(a, 'outlet', 'BCOutflow', 'kmax')
t = C.newPyTree(['Base', a])
bc = Internal.getNodeFromType2(a, 'BC_t')

# On a zone
n = Internal.getBCNodesFromName(a)  # all BCs
test.testO(n, 1)

n = Internal.getBCNodesFromName(a, 'cylinder')  # not found, empty list
test.testO(n, 2)

n = Internal.getBCNodesFromName(a, 'wall1')  # one BC
test.testO(n, 3)

n = Internal.getBCNodesFromName(a, ['wall1', 'wall2'])  # from a list of BC names
test.testO(n, 4)

n = Internal.getBCNodesFromName(a, 'wall*')  # wildcard
test.testO(n, 5)

# On a PyTree
n = Internal.getBCNodesFromName(t)  # all BCs
test.testO(n, 6)

n = Internal.getBCNodesFromName(t, 'cylinder')  # not found, empty list
test.testO(n, 7)

n = Internal.getBCNodesFromName(t, 'wall1')  # one BC
test.testO(n, 8)

n = Internal.getBCNodesFromName(t, 'wall*')  # wildcard
test.testO(n, 9)

# On a BC_t node
n = Internal.getBCNodesFromName(bc)  # all BCs
test.testO(n, 10)

n = Internal.getBCNodesFromName(bc, 'cylinder')  # not found, empty list
test.testO(n, 11)

n = Internal.getBCNodesFromName(bc, 'wall1')  # one BC
test.testO(n, 12)

n = Internal.getBCNodesFromName(bc, 'wall*')  # wildcard
test.testO(n, 13)

# On a list of zones
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((5,0.9,0), (0.01,0.01,1.), (20,20,2)); b[0] = 'cart2'
zones = [a, b]
C._addBC2Zone(zones, 'wall1', 'BCWall', 'imin')
C._addBC2Zone(zones, 'wall2', 'BCWall', 'jmin')
C._addBC2Zone(zones, 'wall3', 'BCWall', 'jmax')
C._addBC2Zone(zones, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(zones, 'farfield', 'BCFarfield', 'kmax')

n = Internal.getBCNodesFromName(zones)  # all BCs
test.testO(n, 14)

n = Internal.getBCNodesFromName(zones, 'cylinder')  # not found, empty list
test.testO(n, 15)

n = Internal.getBCNodesFromName(zones, 'wall3')  # one BC
test.testO(n, 16)

n = Internal.getBCNodesFromName(zones, 'wall*')  # wildcard
test.testO(n, 17)
