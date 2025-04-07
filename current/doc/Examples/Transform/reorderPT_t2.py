# - reorder (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

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
