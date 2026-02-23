# - reorder (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
import numpy

# reorder a Cartesian structured mesh
a = G.cart((0,0,0), (1,1,1), (8,9,20))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
C._initVars(a, 'centers:cellN', 1.)
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imax')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.reorder(t, (2,1,-3))
test.testT(t,1)

# reorder a TRI mesh
a = G.cartTetra((0,0,0), (1,1,1), (8,9,1))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
#C._initVars(a,'centers:cellN',1.)
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.reorder(t, (-1,))
test.testT(t,2)

# reorder a QUAD mesh
a = G.cartHexa((0,0,0), (1,1,1), (8,9,1))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
#C._initVars(a, 'centers:cellN', 1.)
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.reorder(t, (-1,))
test.testT(t,3)

# reorder a TETRA mesh
a = G.cartTetra((0,0,0), (1,1,1), (5,5,5))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
G._getVolumeMap(a)
test.testT(a,4)  # expected results

ec = Internal.getNodeFromName(a, 'ElementConnectivity')[1]
ec = ec.reshape(-1, 4)
ec[:, [1, 2]] = ec[:, [2, 1]]  # flip surface normals by swapping 2 indices
ec = ec.ravel()
G._getVolumeMap(a)  # here negative volumes
test.testT(a,5)

T._reorder(a)
G._getVolumeMap(a)
test.testT(a,4)

# reorder an HEXA mesh
a = G.cartHexa((0,0,0), (1,1,1), (5,5,5))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
G._getVolumeMap(a)
test.testT(a,6)  # expected results

ec = Internal.getNodeFromName(a, 'ElementConnectivity')[1]
ec = ec.reshape(-1, 8)
ec[:, [1, 3]] = ec[:, [3, 1]]  # flip surface normals by swapping 2 indices
ec[:, [5, 7]] = ec[:, [7, 5]]  # and their respective image points
ec = ec.ravel()

G._getVolumeMap(a)  # here negative volumes
test.testT(a,7)

T._reorder(a)
G._getVolumeMap(a)
test.testT(a,6)

# reorder a PYRA mesh
a = G.cartPyra((0,0,0), (1,1,1), (5,5,5))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
G._getVolumeMap(a)
test.testT(a,8)  # expected results

ec = Internal.getNodeFromName(a, 'ElementConnectivity')[1]
ec = ec.reshape(-1, 5)
ec[:, [1, 3]] = ec[:, [3, 1]]  # flip surface normals by swapping 2 indices
ec = ec.ravel()

G._getVolumeMap(a)  # here negative volumes
test.testT(a,9)

T._reorder(a)
G._getVolumeMap(a)
test.testT(a,8)

# reorder a PENTA mesh
a = G.cartPenta((0,0,0), (1,1,1), (5,5,5))
C._initVars(a, '{Density}={CoordinateX}**2+2*{CoordinateY}')
G._getVolumeMap(a)
test.testT(a,10)  # expected results

ec = Internal.getNodeFromName(a, 'ElementConnectivity')[1]
ec = ec.reshape(-1, 6)
ec[:, [1, 2]] = ec[:, [2, 1]]  # flip surface normals by swapping 2 indices
ec[:, [4, 5]] = ec[:, [5, 4]]  # and their respective image points
ec = ec.ravel()

G._getVolumeMap(a)  # here negative volumes
test.testT(a,11)

T._reorder(a)
G._getVolumeMap(a)
test.testT(a,10)
