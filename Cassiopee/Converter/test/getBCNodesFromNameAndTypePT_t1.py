# - getBCNodesFromNameAndType (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (10,10,2))
C._addBC2Zone(a, 'wall1', 'BCWallInviscid', 'jmin')
C._addBC2Zone(a, 'wall2', 'BCWallViscous', 'jmax')
C._addBC2Zone(a, 'farfield', 'BCFarfield', 'kmin')
C._addBC2Zone(a, 'wall3', 'BCWallViscousIsothermal', 'kmax')
C._addBC2Zone(a, 'wall4', 'FamilySpecified:CYLINDER', 'imin')
C._addBC2Zone(a, 'wall5', 'FamilySpecified:CYLINDER', 'imax')
t = C.newPyTree(['Base', a])
C._addFamily2Base(t[2][1], 'CYLINDER', bndType='BCWallViscous')
bc = Internal.getNodeFromType2(a, 'BC_t')

# On a PyTree
n = Internal.getBCNodesFromNameAndType(t)  # all BCs
test.testO(n, 1)

n = Internal.getBCNodesFromNameAndType(t, bndType='BCSymmetryPlane')  # not found, empty list
test.testO(n, 2)

n = Internal.getBCNodesFromNameAndType(t, bndType='BCWallInviscid')  # one BC
test.testO(n, 3)

n = Internal.getBCNodesFromNameAndType(t, bndName=['wall1', 'wall3'])  # from a list of BC names
test.testO(n, 4)

n = Internal.getBCNodesFromNameAndType(t, bndType=['BCFarfield', 'BCWallViscousIsothermal'])  # from a list of BC types
test.testO(n, 5)

n = Internal.getBCNodesFromNameAndType(t, bndName='wall*')  # wildcard on BC names
test.testO(n, 6)

n = Internal.getBCNodesFromNameAndType(t, bndType='BCWall*')  # wildcard on BC types
test.testO(n, 7)

n = Internal.getBCNodesFromNameAndType(t, bndName='wall*', bndType='BCWallViscous')  # combination
test.testO(n, 8)

n = Internal.getBCNodesFromNameAndType(t, bndName='wall*', bndType='BCWallViscous*')  # combination
test.testO(n, 9)

# On a BC_t node that is not FamilySpecified
n = Internal.getBCNodesFromNameAndType(bc)  # all BCs
test.testO(n, 10)

n = Internal.getBCNodesFromNameAndType(bc, bndType='BCWallInviscid')  # one BC
test.testO(n, 11)

n = Internal.getBCNodesFromNameAndType(bc, bndType='BCWall*')  # wildcard
test.testO(n, 12)

# On a list of zones
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((5,0.9,0), (0.01,0.01,1.), (20,20,2)); b[0] = 'cart2'
zones = [a, b]
C._addBC2Zone(zones, 'wall1', 'BCWallInviscid', 'imin')
C._addBC2Zone(zones, 'wall2', 'BCWallViscous', 'jmin')
C._addBC2Zone(zones, 'wall3', 'BCWallViscousIsothermal', 'jmax')
C._addBC2Zone(zones, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(zones, 'farfield', 'BCFarfield', 'kmax')

n = Internal.getBCNodesFromNameAndType(zones)  # all BCs
test.testO(n, 13)

n = Internal.getBCNodesFromNameAndType(zones, bndType='BCSymmetryPlane')  # not found, empty list
test.testO(n, 14)

n = Internal.getBCNodesFromNameAndType(zones, bndName='wall*', bndType='BCWallViscousIsothermal')  # one BC
test.testO(n, 15)

n = Internal.getBCNodesFromNameAndType(zones, bndName='wall*', bndType='BCWallViscous*')  # combination
test.testO(n, 16)
